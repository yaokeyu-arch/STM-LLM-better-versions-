#!/usr/bin/env python3
"""Run background subtraction followed by double-tip analysis."""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
from pathlib import Path

from PIL import Image, ImageDraw


def make_summary(stem: str, output: Path, panels: list[Path], evidence: Path, sliced: Path) -> Path:
    """Build a compact downloadable evidence summary from independently saved products."""
    available = [Image.open(path).convert("RGB") for path in panels if path.exists()]
    if not available:
        raise RuntimeError("No evidence images were produced")
    thumb_size = (420, 320)
    canvas = Image.new("RGB", (1260, 760), "white")
    for index, image in enumerate(available[:6]):
        image.thumbnail(thumb_size)
        x, y = (index % 3) * 420, (index // 3) * 320
        canvas.paste(image, (x, y))
    result = json.loads(evidence.read_text())
    sliced_result = json.loads(sliced.read_text()) if sliced.exists() else {}
    draw = ImageDraw.Draw(canvas)
    text = (
        f"{stem}\nClassification: {result['classification']}  Confidence: {result.get('confidence', 0):.2f}\n"
        f"Sliced evidence: {sliced_result.get('classification', 'not available')}\n"
        f"Reason: {result['reason']}"
    )
    draw.multiline_text((15, 650), text, fill="black", spacing=5)
    summary = output.with_name(f"{stem}_final_evidence_summary.png")
    canvas.save(summary)
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("sxm", type=Path)
    parser.add_argument("-o", "--output", type=Path)
    parser.add_argument("--candidates", type=int, default=10)
    args = parser.parse_args()

    here = Path(__file__).resolve().parent
    output = args.output or Path.cwd() / f"{args.sxm.stem}_background_subtracted.npz"
    evidence = output.with_name(f"{args.sxm.stem}_double_tip_evidence.json")
    acf_evidence = output.with_name(f"{args.sxm.stem}_normalized_acf.npz")
    acf_image = output.with_name(f"{args.sxm.stem}_autocorrelation_marked.png")
    feature_image = output.with_name(f"{args.sxm.stem}_duplicated_features_marked.png")
    original_image = output.with_name(f"{args.sxm.stem}_original.png")
    flattened_image = output.with_name(f"{args.sxm.stem}_flattened.png")
    denoised_image = output.with_name(f"{args.sxm.stem}_noise_reduced.png")
    noise_mask = output.with_name(f"{args.sxm.stem}_noise_mask.png")
    sliced_dir = output.with_name(f"{args.sxm.stem}_sliced")
    subprocess.run(
        [sys.executable, str(here / "stm_background_subtract.py"), str(args.sxm), "-o", str(output)],
        check=True,
    )
    subprocess.run(
        [sys.executable, str(here / "stm_double_tip_analysis.py"), str(args.sxm),
         "--candidates", str(args.candidates), "--json-output", str(evidence),
         "--acf-output", str(acf_evidence), "--acf-image-output", str(acf_image),
         "--feature-image-output", str(feature_image),
         "--original-output", str(original_image),
         "--flattened-output", str(flattened_image),
         "--denoised-output", str(denoised_image),
         "--noise-mask-output", str(noise_mask)],
        check=True,
    )
    subprocess.run(
        [sys.executable, str(here / "Sliced Autocorrection" / "sliced_autocorrelation.py"),
         str(args.sxm), "--output-dir", str(sliced_dir)],
        check=True,
    )
    sliced_summary = sliced_dir / f"{args.sxm.stem}_sliced_summary.json"
    sliced_map = sliced_dir / f"{args.sxm.stem}_sliced_window_map.png"
    summary = make_summary(
        args.sxm.stem, output,
        [original_image, flattened_image, denoised_image, acf_image, feature_image, sliced_map],
        evidence, sliced_summary,
    )
    print(f"Final evidence summary: {summary}")


if __name__ == "__main__":
    main()
