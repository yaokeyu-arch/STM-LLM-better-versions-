#!/usr/bin/env python3
"""Create a non-overwriting STM run layout and print its path."""

from __future__ import annotations

import argparse
import json
import re
from datetime import datetime, timezone
from pathlib import Path


GENERIC_NAMES = {"", ".", "images", "data", "files", "input", "inputs", "untitled"}


def clean_name(value: str) -> str:
    value = re.sub(r"[\\/:*?\"<>|]+", "-", value).strip(" .-")
    return re.sub(r"\s+", " ", value)


def available_numbered_name(root: Path, base: str) -> str:
    if base and base.lower() not in GENERIC_NAMES:
        if not (root / base).exists():
            return base
        number = 2
        while (root / f"{base} - run {number}").exists():
            number += 1
        return f"{base} - run {number}"
    number = 1
    while (root / f"run {number}").exists():
        number += 1
    return f"run {number}"


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source", type=Path, help="folder containing the image group")
    parser.add_argument("runs_root", type=Path, help="workspace folder that stores all runs")
    parser.add_argument("--name", help="explicit run name; defaults to the source-folder name")
    parser.add_argument("--images", nargs="*", default=[], help="image filenames to pre-create")
    args = parser.parse_args()

    args.runs_root.mkdir(parents=True, exist_ok=True)
    requested = clean_name(args.name if args.name is not None else args.source.name)
    run_name = available_numbered_name(args.runs_root, requested)
    run_dir = args.runs_root / run_name
    overall_dir = run_dir / "overall-results"
    images_dir = run_dir / "images"
    overall_dir.mkdir(parents=True)
    images_dir.mkdir()

    image_entries = []
    for filename in args.images:
        image_name = Path(filename).stem
        image_dir = images_dir / image_name
        image_dir.mkdir(exist_ok=True)
        image_entries.append({"input": filename, "output_folder": str(image_dir.relative_to(run_dir))})

    manifest = {
        "protocol": "saving-protocol",
        "run_name": run_name,
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "status": "created",
        "source_group": str(args.source.resolve()),
        "overall_results_folder": "overall-results",
        "images_folder": "images",
        "image_count": len(image_entries),
        "images": image_entries,
    }
    (run_dir / "run_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(run_dir.resolve())


if __name__ == "__main__":
    main()
