#!/usr/bin/env python3
"""Background-subtract Nanonis SXM Z data before double-tip analysis.

The operation removes a fitted plane and isolated scan-line offsets, but does
not flatten away local morphology. Both forward and backward Z scans are saved
in physical height units to a compressed NPZ file.

Usage:
    python stm_background_subtract.py scan.sxm
    python stm_double_tip_analysis.py scan_background_subtracted.npz
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np

from stm_double_tip_analysis import background_subtract, crop_common_finite, load_sxm_z


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("sxm", type=Path, help="input Nanonis .sxm file")
    parser.add_argument("-o", "--output", type=Path, help="output .npz path")
    args = parser.parse_args()

    forward, backward, pixel_x, pixel_y = load_sxm_z(args.sxm)
    original_shape = forward.shape
    forward, backward, crop = crop_common_finite(forward, backward)
    forward = background_subtract(forward)
    backward = background_subtract(backward)
    output = args.output or args.sxm.with_name(
        f"{args.sxm.stem}_background_subtracted.npz"
    )
    output.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        output,
        z_forward=forward,
        z_backward=backward,
        pixel_x_m=pixel_x,
        pixel_y_m=pixel_y,
        source=str(args.sxm),
        preprocessing="least-squares plane plus robust isolated row-offset correction",
        finite_crop=np.asarray(crop),
        original_shape=np.asarray(original_shape),
    )
    print(f"Saved: {output}")
    print(f"Forward RMS height: {np.std(forward) * 1e12:.3f} pm")
    print(f"Backward RMS height: {np.std(backward) * 1e12:.3f} pm")
    if forward.shape != original_shape:
        print(f"Finite-data crop: x={crop[0]}:{crop[2]}, y={crop[1]}:{crop[3]}")


if __name__ == "__main__":
    main()
