#!/usr/bin/env python3
"""Run the packaged STM analysis skill with a reproducible saving archive.

Examples:
  python3 main.py scan.sxm --run-name "scan verification"
  python3 main.py /path/to/scans --run-name "sample group"
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import importlib.metadata
import json
import os
import platform
import shutil
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path


SCRIPT_DIR = Path(__file__).resolve().parent
SKILL_DIR = SCRIPT_DIR.parent
WORKSPACE = SKILL_DIR.parents[1]
SAVING_SKILL = WORKSPACE / "skills" / "saving-protocol"
DEFAULT_RUNS = WORKSPACE / "runs"


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def append_command(path: Path, command: list[str], cwd: Path, started: str,
                   exit_code: int, stdout: str, stderr: str) -> None:
    record = {
        "started_utc": started,
        "finished_utc": utc_now(),
        "cwd": str(cwd),
        "command": command,
        "exit_code": exit_code,
        "stdout_tail": stdout[-4000:],
        "stderr_tail": stderr[-4000:],
    }
    with path.open("a", encoding="utf-8") as stream:
        stream.write(json.dumps(record) + "\n")


def run_logged(command: list[str], cwd: Path, commands_log: Path) -> subprocess.CompletedProcess:
    started = utc_now()
    result = subprocess.run(command, cwd=cwd, text=True, capture_output=True)
    append_command(commands_log, command, cwd, started, result.returncode,
                   result.stdout, result.stderr)
    return result


def create_layout(source_group: Path, runs_root: Path, run_name: str | None,
                  images: list[Path]) -> Path:
    helper = SAVING_SKILL / "scripts" / "create_run_layout.py"
    command = [sys.executable, str(helper), str(source_group), str(runs_root)]
    if run_name:
        command += ["--name", run_name]
    command += ["--images", *[image.name for image in images]]
    result = subprocess.run(command, text=True, capture_output=True, check=True)
    return Path(result.stdout.strip().splitlines()[-1]).resolve()


def archive_materials(run_dir: Path) -> None:
    code_dir = run_dir / "materials" / "code"
    skill_dir = run_dir / "materials" / "skills"
    env_dir = run_dir / "materials" / "environment"
    code_dir.mkdir(parents=True, exist_ok=True)
    skill_dir.mkdir(parents=True, exist_ok=True)
    env_dir.mkdir(parents=True, exist_ok=True)

    for relative in (
        "main.py", "stm_background_subtract.py", "stm_double_tip_analysis.py",
        "run_stm_double_tip_pipeline.py", "batch_analyze_folder.py",
        "Sliced Autocorrection/sliced_autocorrelation.py",
    ):
        source = SCRIPT_DIR / relative
        if source.exists():
            destination = code_dir / relative
            destination.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(source, destination)
    shutil.copy2(SKILL_DIR / "SKILL.md", skill_dir / "stm-double-tip-analysis.SKILL.md")
    shutil.copy2(SAVING_SKILL / "SKILL.md", skill_dir / "saving-protocol.SKILL.md")

    packages = {}
    for name in ("numpy", "scipy", "Pillow"):
        try:
            packages[name] = importlib.metadata.version(name)
        except importlib.metadata.PackageNotFoundError:
            packages[name] = None
    environment = {
        "recorded_utc": utc_now(),
        "python": sys.version,
        "executable": sys.executable,
        "platform": platform.platform(),
        "packages": packages,
    }
    (env_dir / "runtime.json").write_text(json.dumps(environment, indent=2) + "\n")


def archive_sources(run_dir: Path, images: list[Path]) -> list[dict]:
    archived = []
    for source in images:
        image_dir = run_dir / "images" / source.stem
        source_dir = image_dir / "source"
        source_dir.mkdir(parents=True, exist_ok=True)
        destination = source_dir / source.name
        shutil.copy2(source, destination)
        source_hash = sha256(source)
        copied_hash = sha256(destination)
        if source_hash != copied_hash:
            raise RuntimeError(f"checksum mismatch after copying {source}")
        archived.append({
            "input": source.name,
            "external_source": str(source),
            "archived_source": str(destination.relative_to(run_dir)),
            "bytes": source.stat().st_size,
            "sha256": source_hash,
        })
    return archived


def run_single(source: Path, run_dir: Path, commands_log: Path,
               candidates: int) -> tuple[bool, str]:
    image_dir = run_dir / "images" / source.stem
    image_dir.mkdir(parents=True, exist_ok=True)
    output = image_dir / f"{source.stem}_background_subtracted.npz"
    command = [
        sys.executable, str(SCRIPT_DIR / "run_stm_double_tip_pipeline.py"),
        str(source), "-o", str(output), "--candidates", str(candidates),
    ]
    result = run_logged(command, WORKSPACE, commands_log)
    (image_dir / "pipeline_stdout.log").write_text(result.stdout)
    (image_dir / "pipeline_stderr.log").write_text(result.stderr)
    return result.returncode == 0, result.stderr[-1000:]


def build_single_summary(run_dir: Path, source: Path) -> None:
    image_dir = run_dir / "images" / source.stem
    evidence_path = image_dir / f"{source.stem}_double_tip_evidence.json"
    sliced_path = image_dir / f"{source.stem}_sliced" / f"{source.stem}_sliced_summary.json"
    evidence = json.loads(evidence_path.read_text())
    sliced = json.loads(sliced_path.read_text())
    row = {
        "file": source.name,
        "whole_classification": evidence["classification"],
        "confidence": evidence["confidence"],
        "whole_reason": evidence["reason"],
        "sliced_classification": sliced["classification"],
        "max_sliced_score": max(sliced.get("scores", [0])),
    }
    overall = run_dir / "overall-results"
    overall.mkdir(exist_ok=True)
    (overall / "combined_analysis.json").write_text(json.dumps([row], indent=2) + "\n")
    with (overall / "combined_analysis.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=row.keys())
        writer.writeheader()
        writer.writerow(row)


def write_inventory(run_dir: Path) -> None:
    entries = []
    for path in sorted(run_dir.rglob("*")):
        if not path.is_file() or path.name == "file_inventory.json":
            continue
        relative = path.relative_to(run_dir)
        role = "analysis evidence"
        if relative.parts[0] == "materials":
            role = "reproducibility material"
        elif "source" in relative.parts:
            role = "archived original input"
        elif relative.parts[0] == "overall-results":
            role = "group-level result"
        entries.append({
            "path": str(relative), "bytes": path.stat().st_size,
            "sha256": sha256(path), "role": role,
        })
    payload = {"created_utc": utc_now(), "file_count": len(entries), "files": entries}
    (run_dir / "file_inventory.json").write_text(json.dumps(payload, indent=2) + "\n")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source", type=Path, help="one SXM file or a folder of SXM files")
    parser.add_argument("--run-name")
    parser.add_argument("--runs-root", type=Path, default=DEFAULT_RUNS)
    parser.add_argument("--candidates", type=int, default=10)
    args = parser.parse_args()

    source = args.source.expanduser().resolve()
    images = sorted(source.glob("*.sxm")) if source.is_dir() else [source]
    images = [path for path in images if path.is_file() and path.suffix.lower() == ".sxm"]
    if not images:
        parser.error("no SXM files found")
    source_group = source if source.is_dir() else source.parent
    run_dir = create_layout(source_group, args.runs_root.resolve(), args.run_name, images)
    materials = run_dir / "materials"
    materials.mkdir(exist_ok=True)
    commands_log = materials / "commands.jsonl"
    archive_materials(run_dir)
    archived_sources = archive_sources(run_dir, images)

    successes, failures = 0, []
    if source.is_dir():
        command = [sys.executable, str(SCRIPT_DIR / "batch_analyze_folder.py"),
                   str(source), str(run_dir)]
        result = run_logged(command, WORKSPACE, commands_log)
        successes = len(images) if result.returncode == 0 else 0
        if result.returncode:
            failures.append({"source": str(source), "error": result.stderr[-1000:]})
    else:
        ok, error = run_single(source, run_dir, commands_log, args.candidates)
        successes = int(ok)
        if ok:
            build_single_summary(run_dir, source)
        else:
            failures.append({"source": str(source), "error": error})

    manifest_path = run_dir / "run_manifest.json"
    manifest = json.loads(manifest_path.read_text())
    manifest.update({
        "status": "complete" if not failures else "complete_with_failures",
        "completed_utc": utc_now(),
        "successful_images": successes,
        "failed_images": len(failures),
        "failures": failures,
        "archived_sources": archived_sources,
        "materials": {
            "code": "materials/code", "skills": "materials/skills",
            "environment": "materials/environment/runtime.json",
            "commands": "materials/commands.jsonl",
        },
        "overall_results": [
            "overall-results/combined_analysis.json",
            "overall-results/combined_analysis.csv",
        ],
        "inventory": "file_inventory.json",
    })
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n")
    write_inventory(run_dir)
    print(run_dir)
    return 0 if not failures else 1


if __name__ == "__main__":
    raise SystemExit(main())
