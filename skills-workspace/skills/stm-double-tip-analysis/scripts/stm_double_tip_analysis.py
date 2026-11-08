#!/usr/bin/env python3
"""Evidence-based double/multiple-tip screening for Nanonis SXM topographs.

Usage:
    python stm_double_tip_analysis.py scan.sxm

Requires NumPy and SciPy. Coordinates are reported with +x right and +y down.
The script uses symmetric autocorrelation peaks only to generate candidates;
classification also requires repeated overlap of irregular morphology.
"""

from __future__ import annotations

import argparse
import json
import math
import re
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from PIL import Image, ImageDraw
from scipy import ndimage


@dataclass
class Candidate:
    dx: int
    dy: int
    acf_forward: float
    acf_backward: float
    prominence_z: float
    overlap_fraction: float
    corr_forward: float
    corr_backward: float
    f1_forward: float
    f1_backward: float
    consistent_tiles_forward: int
    consistent_tiles_backward: int
    feature_tiles_forward: int
    feature_tiles_backward: int
    periodicity_penalty: float
    translation_corr_forward: float
    translation_corr_backward: float
    translation_f1_forward: float
    translation_f1_backward: float
    scale_forward: float
    scale_backward: float
    rotation_forward_deg: float
    rotation_backward_deg: float
    mirrored_forward: bool
    mirrored_backward: bool


def _field(header: str, name: str) -> str:
    match = re.search(rf":{re.escape(name)}:\s*\n([^\n]+)", header)
    if not match:
        raise ValueError(f"Missing SXM header field: {name}")
    return match.group(1).strip()


def load_sxm_z(path: Path) -> tuple[np.ndarray, np.ndarray, float, float]:
    raw = path.read_bytes()
    marker = b":SCANIT_END:"
    marker_at = raw.find(marker)
    if marker_at < 0:
        raise ValueError("Not a recognized Nanonis SXM file")

    header = raw[: marker_at + len(marker)].decode("latin1", "replace")
    nx, ny = map(int, _field(header, "SCAN_PIXELS").split()[:2])
    range_x, range_y = map(float, _field(header, "SCAN_RANGE").split()[:2])

    info_match = re.search(r":DATA_INFO:\s*\n(.*?)(?=\n:SCANIT_END:)", header, re.S)
    if not info_match:
        raise ValueError("Missing DATA_INFO table")
    rows = [line.strip() for line in info_match.group(1).splitlines() if line.strip()]
    channels: list[tuple[str, str]] = []
    for row in rows[1:]:
        columns = row.split("\t")
        if len(columns) >= 4:
            channels.append((columns[1].strip(), columns[3].strip().lower()))

    z_index = next((i for i, (name, _) in enumerate(channels) if name == "Z"), None)
    if z_index is None:
        raise ValueError("No Z channel found")

    frame_counts = [2 if direction == "both" else 1 for _, direction in channels]
    frame_count = sum(frame_counts)
    data_bytes = frame_count * nx * ny * 4
    data_offset = len(raw) - data_bytes
    if data_offset <= marker_at:
        raise ValueError("SXM data size is inconsistent with its header")

    frames = np.frombuffer(raw, dtype=">f4", offset=data_offset).reshape(frame_count, ny, nx)
    first = sum(frame_counts[:z_index])
    forward = frames[first].astype(float)
    if frame_counts[z_index] == 2:
        # Nanonis stores the backward trace in acquisition order.
        backward = frames[first + 1][:, ::-1].astype(float)
    else:
        backward = forward.copy()
    return forward, backward, range_x / nx, range_y / ny


def load_sxm_metadata(path: Path) -> dict:
    """Read acquisition metadata without using it to influence classification."""
    if path.suffix.lower() == ".npz":
        with np.load(path) as data:
            if "source" in data.files:
                source = Path(str(data["source"].item()))
                if source.exists():
                    return load_sxm_metadata(source)
            metadata = {}
            if "original_shape" in data.files:
                ny, nx = map(int, data["original_shape"])
                metadata.update({"pixels_x": nx, "pixels_y": ny})
            if "pixel_x_m" in data.files and "pixel_y_m" in data.files:
                metadata["pixel_size_m"] = {
                    "x": float(data["pixel_x_m"]), "y": float(data["pixel_y_m"])
                }
            return metadata
    if path.suffix.lower() != ".sxm":
        return {}
    raw = path.read_bytes()
    marker = b":SCANIT_END:"
    marker_at = raw.find(marker)
    if marker_at < 0:
        return {}
    header = raw[: marker_at + len(marker)].decode("latin1", "replace")
    nx, ny = map(int, _field(header, "SCAN_PIXELS").split()[:2])
    range_x, range_y = map(float, _field(header, "SCAN_RANGE").split()[:2])
    try:
        scan_direction = _field(header, "SCAN_DIR").lower()
    except ValueError:
        scan_direction = "unknown"
    return {
        "pixels_x": nx,
        "pixels_y": ny,
        "scan_range_m": {"x": range_x, "y": range_y},
        "pixel_size_m": {"x": range_x / nx, "y": range_y / ny},
        "slow_scan_direction": scan_direction,
        "forward_backward_available": True,
    }


def background_subtract(z: np.ndarray) -> np.ndarray:
    """Remove a plane and isolated row offsets while preserving local features."""
    h, w = z.shape
    y, x = np.mgrid[:h, :w]
    design = np.c_[np.ones(z.size), x.ravel(), y.ravel()]
    plane = design @ np.linalg.lstsq(design, z.ravel(), rcond=None)[0]
    residual = z - plane.reshape(h, w)

    row_medians = np.median(residual, axis=1)
    trend = ndimage.median_filter(row_medians, size=max(5, h // 15 | 1))
    residual -= (row_medians - trend)[:, None]
    return residual


NOISE_PARAMETERS = {
    "isolated_glitch_threshold_mad_sigma": 10.0,
    "maximum_glitch_component_pixels": 2,
    "replacement": "3x3 local median only inside confirmed isolated-glitch mask",
    "general_smoothing": "none",
}


def conservative_noise_reduce(flattened: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Remove only isolated extreme glitches; preserve small and repeated morphology."""
    local_median = ndimage.median_filter(flattened, size=3, mode="nearest")
    deviation = flattened - local_median
    center = float(np.median(deviation))
    robust_sigma = 1.4826 * float(np.median(np.abs(deviation - center))) + 1e-12
    extreme = np.abs(deviation - center) > (
        NOISE_PARAMETERS["isolated_glitch_threshold_mad_sigma"] * robust_sigma
    )
    labels, count = ndimage.label(extreme, structure=np.ones((3, 3), dtype=bool))
    sizes = np.bincount(labels.ravel(), minlength=count + 1)
    isolated_labels = np.flatnonzero(
        (sizes <= NOISE_PARAMETERS["maximum_glitch_component_pixels"]) & (np.arange(len(sizes)) > 0)
    )
    noise_mask = np.isin(labels, isolated_labels)
    reduced = flattened.copy()
    reduced[noise_mask] = local_median[noise_mask]
    return reduced, noise_mask


def preprocess_with_evidence(z: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Return flattened, noise-reduced, noise-mask, and normalized analysis maps."""
    flattened = background_subtract(z)
    reduced, noise_mask = conservative_noise_reduce(flattened)
    scale = float(np.std(reduced))
    analysis = (reduced - float(np.mean(reduced))) / max(scale, 1e-12)
    return flattened, reduced, noise_mask, analysis


def preprocess(z: np.ndarray) -> np.ndarray:
    """Compatibility wrapper returning the normalized noise-reduced analysis map."""
    return preprocess_with_evidence(z)[3]


def save_noise_removed_image(forward: np.ndarray, backward: np.ndarray, path: Path) -> None:
    """Save a robustly scaled view after conservative morphology-preserving denoising."""
    denoised = (forward + backward) / 2
    low, high = np.quantile(denoised[np.isfinite(denoised)], [0.01, 0.99])
    scaled = np.clip((denoised - low) / max(high - low, 1e-12), 0, 1)
    pixels = np.round(255 * scaled).astype(np.uint8)
    path.parent.mkdir(parents=True, exist_ok=True)
    Image.fromarray(pixels, mode="L").save(path)


def save_noise_mask(mask: np.ndarray, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    Image.fromarray(mask.astype(np.uint8) * 255, mode="L").save(path)


def grayscale_image(image: np.ndarray) -> Image.Image:
    finite = image[np.isfinite(image)]
    low, high = np.quantile(finite, [0.01, 0.99])
    scaled = np.clip((image - low) / max(high - low, 1e-12), 0, 1)
    return Image.fromarray(np.round(255 * scaled).astype(np.uint8), mode="L")


def orbital_lobe_points(
    image: np.ndarray, support_quantile: float = 0.80
) -> tuple[np.ndarray, np.ndarray, np.ndarray, list[tuple[int, int]]]:
    """Return individual bright-lobe centers, including weaker touching lobes.

    A lightly smoothed map preserves a weaker shoulder that can disappear in
    the broader morphology image. The support mask remains height selective,
    and a robust intensity floor prevents noise maxima from becoming lobes.
    """
    fine = ndimage.gaussian_filter(image, 0.35)
    smooth = ndimage.gaussian_filter(image, 0.70)
    mask = smooth >= np.quantile(smooth, support_quantile)
    center = float(np.median(fine))
    robust_sigma = 1.4826 * float(np.median(np.abs(fine - center))) + 1e-12
    supported = ndimage.binary_dilation(mask, structure=np.ones((3, 3), dtype=bool))
    maxima = (
        (fine == ndimage.maximum_filter(fine, size=3))
        & supported
        & (fine >= center + 0.60 * robust_sigma)
    )
    maximum_labels, maximum_count = ndimage.label(maxima)
    points: list[tuple[int, int]] = []
    for maximum_id in range(1, maximum_count + 1):
        plateau = maximum_labels == maximum_id
        values = np.where(plateau, fine, -np.inf)
        points.append(tuple(map(int, np.unravel_index(np.argmax(values), values.shape))))

    # A touching doublet can have only one discrete 3x3 maximum even though
    # its bright core is clearly elongated into two orbital lobes.  Deblend
    # compact elongated cores along their principal axis and add the two
    # weighted half-centres as lobe candidates.  These are candidates only;
    # the later one-to-one displacement and morphology tests still decide
    # whether they form accepted evidence.
    core = fine >= max(float(np.quantile(fine, 0.86)), center + 0.80 * robust_sigma)
    core_labels, core_count = ndimage.label(core, structure=np.ones((3, 3), dtype=bool))
    for core_id in range(1, core_count + 1):
        yy, xx = np.where(core_labels == core_id)
        area = len(xx)
        if not 8 <= area <= 250:
            continue
        weights = np.maximum(fine[yy, xx] - center, 1e-12)
        coordinates = np.column_stack((xx, yy)).astype(float)
        centroid = np.average(coordinates, axis=0, weights=weights)
        centered = coordinates - centroid
        covariance = (centered * weights[:, None]).T @ centered / max(float(weights.sum()), 1e-12)
        eigenvalues, eigenvectors = np.linalg.eigh(covariance)
        if eigenvalues[0] <= 1e-12 or eigenvalues[1] / eigenvalues[0] < 2.0:
            continue
        axis = eigenvectors[:, 1]
        projection = centered @ axis
        if float(projection.max() - projection.min()) < 4.5:
            continue
        left = projection < 0
        right = ~left
        if left.sum() < 3 or right.sum() < 3:
            continue
        centres: list[tuple[int, int]] = []
        for half in (left, right):
            half_xy = np.average(coordinates[half], axis=0, weights=weights[half])
            hx = int(np.clip(round(float(half_xy[0])), 0, image.shape[1] - 1))
            hy = int(np.clip(round(float(half_xy[1])), 0, image.shape[0] - 1))
            centres.append((hy, hx))
        separation = math.hypot(centres[1][1] - centres[0][1], centres[1][0] - centres[0][0])
        if not 2.0 <= separation <= 14.0:
            continue
        for candidate in centres:
            if not any(math.hypot(candidate[1] - px, candidate[0] - py) <= 1.5 for py, px in points):
                points.append(candidate)
    labels, _ = ndimage.label(mask, structure=np.ones((3, 3), dtype=bool))
    return smooth, mask, labels, points


def repeated_close_lobe_vectors(image: np.ndarray) -> list[tuple[int, int, int]]:
    """Find short vectors repeated by resolved or touching bright lobes.

    This supplements ACF maxima when a real doublet lies on the shoulder of
    the central ACF peak and therefore is not itself a resolved ACF maximum.
    Only vectors repeated in several separate bright components are returned.
    """
    smooth, mask, components, all_points = orbital_lobe_points(image)
    component_count = int(components.max())
    maxima_points: dict[int, list[tuple[int, int]]] = {}
    for ay, ax in all_points:
        component_id = int(components[ay, ax])
        if not component_id:
            nearby = components[max(0, ay-1):ay+2, max(0, ax-1):ax+2]
            nonzero = nearby[nearby > 0]
            component_id = int(nonzero[0]) if len(nonzero) else 0
        if component_id:
            maxima_points.setdefault(component_id, []).append((ay, ax))
    counts: dict[tuple[int, int], set[int]] = {}
    for component_id in range(1, component_count + 1):
        points = maxima_points.get(component_id, [])
        if len(points) < 2:
            continue
        for index, (ay, ax) in enumerate(points):
            for by, bx in points[index + 1:]:
                dx, dy = bx - ax, by - ay
                distance = math.hypot(dx, dy)
                if not 2.0 <= distance <= 14.0:
                    continue
                height_similarity = min(smooth[ay, ax], smooth[by, bx]) / max(
                    smooth[ay, ax], smooth[by, bx], 1e-12
                )
                if height_similarity < 0.45:
                    continue
                if dy < 0 or (dy == 0 and dx < 0):
                    dx, dy = -dx, -dy
                counts.setdefault((int(dx), int(dy)), set()).add(component_id)
    return [
        (dx, dy, len(component_ids))
        for (dx, dy), component_ids in counts.items()
        if len(component_ids) >= 3
    ]


def feature_pairs(image: np.ndarray, dx: int, dy: int, tolerance: float = 3.0) -> list[dict]:
    """Locate reciprocal one-to-one lobe pairs at the candidate displacement.

    Each local maximum is an individual lobe, including maxima inside one
    touching bright component. A lobe may participate in at most one pair
    within this vector family. Final evidence must additionally arbitrate
    competing assignments across vector families; a short nearby match must
    not supersede a better long-range orbital match merely because it is
    closer or has more similar absolute brightness.
    """
    smooth, mask, labels, points = orbital_lobe_points(image)
    effective_tolerance = min(tolerance, max(1.25, 0.20 * math.hypot(dx, dy)))

    radius = 6
    proposed: list[tuple[float, float, float, tuple[int, int], tuple[int, int]]] = []
    for ay, ax in points:
        if not (radius <= ax < image.shape[1] - radius and radius <= ay < image.shape[0] - radius):
            continue
        patch_a = image[ay-radius:ay+radius+1, ax-radius:ax+radius+1]
        target_x, target_y = ax + dx, ay + dy
        for by, bx in points:
            if (ay, ax) == (by, bx):
                continue
            vector_error = math.hypot(bx - target_x, by - target_y)
            if vector_error > effective_tolerance:
                continue
            if not (radius <= bx < image.shape[1] - radius and radius <= by < image.shape[0] - radius):
                continue
            patch_b = image[by-radius:by+radius+1, bx-radius:bx+radius+1]
            translation = float(np.corrcoef(patch_a.ravel(), patch_b.ravel())[0, 1])
            if not np.isfinite(translation):
                continue
            # Displacement agreement dominates; morphology breaks ties between
            # multiple lobes near the expected partner location.
            assignment_score = 0.70 * max(translation, 0.0) + 0.30 * (1.0 - vector_error / effective_tolerance)
            proposed.append((assignment_score, vector_error, translation, (ay, ax), (by, bx)))

    # Greedy maximum-score set packing enforces reciprocal one-to-one use of
    # lobes on both sides.  A point cannot be reused as A or B in another pair.
    proposed.sort(key=lambda item: (item[0], -item[1]), reverse=True)
    selected: list[tuple[float, tuple[int, int], tuple[int, int]]] = []
    used: set[tuple[int, int]] = set()
    for _, _, translation, point_a, point_b in proposed:
        if point_a in used or point_b in used:
            continue
        selected.append((translation, point_a, point_b))
        used.update((point_a, point_b))

    pairs: list[dict] = []
    for pair_index, (translation, (ay, ax), (by, bx)) in enumerate(selected, 1):
        patch_a = image[ay-radius:ay+radius+1, ax-radius:ax+radius+1]
        patch_b = image[by-radius:by+radius+1, bx-radius:bx+radius+1]
        best = (translation, 1.0, 0.0, False)
        for scale in (0.90, 0.95, 1.0, 1.05, 1.10):
            for rotation in (-6.0, -3.0, 0.0, 3.0, 6.0):
                for mirrored in (False, True):
                    transformed = transformed_patch(patch_b, scale, rotation, mirrored)
                    score = float(np.corrcoef(patch_a.ravel(), transformed.ravel())[0, 1])
                    penalty = 0.10 * abs(scale - 1.0) / 0.10 + 0.08 * abs(rotation) / 6.0 + 0.04 * mirrored
                    if score - penalty > best[0]:
                        best = (score, scale, rotation, mirrored)
        def component_at(y: int, x: int) -> int:
            component = int(labels[y, x])
            if component:
                return component
            nearby = labels[max(0, y-1):y+2, max(0, x-1):x+2]
            nonzero = nearby[nearby > 0]
            return int(nonzero[0]) if len(nonzero) else 0

        component_a, component_b = component_at(ay, ax), component_at(by, bx)
        area_a = int(np.count_nonzero(labels == component_a)) if component_a else 1
        area_b = int(np.count_nonzero(labels == component_b)) if component_b else 1
        shape_similarity = min(area_a, area_b) / max(area_a, area_b, 1)
        height_similarity = min(abs(float(image[ay, ax])), abs(float(image[by, bx]))) / max(
            abs(float(image[ay, ax])), abs(float(image[by, bx])), 1e-12
        )
        deformation = float(max(0.0, 1.0 - shape_similarity))
        same_bright_component = bool(component_a != 0 and component_a == component_b)
        valid_pair = bool(
            translation >= 0.30
            and best[0] >= 0.35
            and height_similarity >= 0.45
            and (shape_similarity >= 0.50 or same_bright_component)
            and (deformation <= 0.50 or same_bright_component)
        )
        if not valid_pair:
            continue
        confidence = float(np.clip(0.45 * max(translation, 0) + 0.25 * max(best[0], 0) + 0.20 * shape_similarity + 0.10 * height_similarity, 0, 1))
        pairs.append({
            "pair_id": "",
            "feature_a_xy": [ax, ay], "feature_b_xy": [bx, by],
            "dx": bx - ax, "dy": by - ay,
            "height_a": float(image[ay, ax]), "height_b": float(image[by, bx]),
            "size_a_pixels": area_a, "size_b_pixels": area_b,
            "shape_similarity": float(shape_similarity),
            "height_similarity": float(height_similarity),
            "translation_only_score": translation,
            "transformation_tolerant_score": float(best[0]),
            "scale_ratio": float(best[1]),
            "rotation_degrees": float(best[2]),
            "mirrored": bool(best[3]),
            "deformation_amount": deformation,
            "touch_or_overlap": same_bright_component,
            "one_to_one_reciprocal": True,
            "pair_validation": "accepted: morphology, height, deformation, and displacement checks passed",
            "vector_tolerance_pixels": float(effective_tolerance),
            "confidence": confidence,
        })
    pairs.sort(key=lambda item: item["confidence"], reverse=True)
    for pair_index, pair in enumerate(pairs, 1):
        pair["pair_id"] = f"P{pair_index:02d}"
    return pairs[:50]


def save_acf_image(acf: np.ndarray, candidates_found: list[Candidate], path: Path) -> None:
    canvas = grayscale_image(acf).convert("RGB")
    draw = ImageDraw.Draw(canvas)
    cy, cx = (acf.shape[0] - 1) // 2, (acf.shape[1] - 1) // 2
    draw.ellipse((cx-3, cy-3, cx+3, cy+3), outline=(255, 220, 0), width=2)
    for index, item in enumerate(candidates_found):
        color = (255, 40, 40) if index == 0 else (0, 220, 255)
        for px, py in ((cx + item.dx, cy + item.dy), (cx - item.dx, cy - item.dy)):
            draw.ellipse((px-4, py-4, px+4, py+4), outline=color, width=2)
    path.parent.mkdir(parents=True, exist_ok=True)
    canvas.save(path)


def save_feature_pair_image(image: np.ndarray, pairs: list[dict], path: Path) -> None:
    canvas = grayscale_image(image).convert("RGB")
    draw = ImageDraw.Draw(canvas)
    colors = [(255, 70, 70), (0, 230, 255), (255, 220, 0), (120, 255, 90), (255, 100, 230)]
    for index, pair in enumerate(pairs[:25]):
        ax, ay = pair["feature_a_xy"]
        bx, by = pair["feature_b_xy"]
        color = colors[index % len(colors)]
        label = pair.get("pair_id", f"P{index + 1:02d}")
        draw.ellipse((ax-3, ay-3, ax+3, ay+3), outline=color, width=2)
        draw.ellipse((bx-3, by-3, bx+3, by+3), outline=color, width=2)
        draw.line((ax, ay, bx, by), fill=color, width=1)
        draw.text((ax + 4, ay - 8), f"{label}A", fill=color)
        draw.text((bx + 4, by - 8), f"{label}B", fill=color)
    path.parent.mkdir(parents=True, exist_ok=True)
    canvas.save(path)


def acf_peak_geometry(acf: np.ndarray, dx: int, dy: int) -> dict:
    cy, cx = (acf.shape[0] - 1) // 2, (acf.shape[1] - 1) // 2
    py, px = cy + dy, cx + dx
    radius = 7
    local = acf[max(0, py-radius):py+radius+1, max(0, px-radius):px+radius+1]
    baseline = float(np.median(local))
    threshold = baseline + 0.5 * (float(acf[py, px]) - baseline)
    mask = local >= threshold
    labels, _ = ndimage.label(mask)
    ly, lx = min(radius, local.shape[0]-1), min(radius, local.shape[1]-1)
    region = labels == labels[ly, lx] if labels[ly, lx] else mask
    coords = np.argwhere(region)
    if len(coords) >= 2:
        covariance = np.cov(coords.T)
        eigenvalues = np.sort(np.linalg.eigvalsh(covariance))[::-1]
        eccentricity = float(np.sqrt(max(0.0, 1.0 - eigenvalues[-1] / max(eigenvalues[0], 1e-12))))
    else:
        eccentricity = 0.0
    positive = float(acf[py, px])
    negative = float(acf[cy-dy, cx-dx])
    return {
        "radial_displacement_pixels": float(math.hypot(dx, dy)),
        "direction_degrees": float(math.degrees(math.atan2(dy, dx))),
        "symmetric_pair_error": float(abs(positive-negative) / max(abs(positive), abs(negative), 1e-12)),
        "peak_size_pixels": int(np.count_nonzero(region)),
        "peak_eccentricity": eccentricity,
        "peak_shape": "elongated/orbital" if eccentricity >= 0.65 else "compact/orbital",
    }


def load_input(path: Path) -> tuple[np.ndarray, np.ndarray, float, float]:
    """Load either a raw SXM file or background-subtracted NPZ output."""
    if path.suffix.lower() == ".sxm":
        return load_sxm_z(path)
    if path.suffix.lower() == ".npz":
        with np.load(path) as data:
            required = {"z_forward", "z_backward", "pixel_x_m", "pixel_y_m"}
            missing = required.difference(data.files)
            if missing:
                raise ValueError(f"NPZ is missing fields: {', '.join(sorted(missing))}")
            return (
                data["z_forward"].astype(float),
                data["z_backward"].astype(float),
                float(data["pixel_x_m"]),
                float(data["pixel_y_m"]),
            )
    raise ValueError("Input must be a Nanonis .sxm file or background-subtracted .npz file")


def crop_common_finite(
    forward: np.ndarray, backward: np.ndarray
) -> tuple[np.ndarray, np.ndarray, tuple[int, int, int, int]]:
    """Keep the largest common fully finite rectangle from an interrupted scan."""
    if forward.shape != backward.shape:
        raise ValueError("Forward and backward Z maps have different shapes")
    valid = np.isfinite(forward) & np.isfinite(backward)
    good_rows = np.all(valid, axis=1)

    def longest_run(mask: np.ndarray) -> tuple[int, int]:
        best_start = best_end = start = 0
        inside = False
        for index, value in enumerate(np.r_[mask, False]):
            if value and not inside:
                start, inside = index, True
            elif not value and inside:
                if index - start > best_end - best_start:
                    best_start, best_end = start, index
                inside = False
        return best_start, best_end

    y0, y1 = longest_run(good_rows)
    if y1 <= y0:
        raise ValueError("No complete rows shared by the forward and backward scans")
    good_cols = np.all(valid[y0:y1], axis=0)
    x0, x1 = longest_run(good_cols)
    retained = (y1 - y0) * (x1 - x0) / valid.size
    if x1 <= x0 or retained < 0.20:
        raise ValueError("Too little common finite Z data for reliable analysis")
    return forward[y0:y1, x0:x1], backward[y0:y1, x0:x1], (x0, y0, x1, y1)


def normalized_acf(image: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    h, w = image.shape
    shape = (2 * h - 1, 2 * w - 1)
    axes = (0, 1)
    spectrum = np.fft.rfftn(image, shape, axes=axes)
    numerator = np.fft.fftshift(
        np.fft.irfftn(spectrum * np.conj(spectrum), shape, axes=axes)
    )
    mask_spectrum = np.fft.rfftn(np.ones_like(image), shape, axes=axes)
    overlap = np.fft.fftshift(
        np.fft.irfftn(mask_spectrum * np.conj(mask_spectrum), shape, axes=axes)
    )
    acf = numerator / np.maximum(overlap, 1)
    acf /= acf[h - 1, w - 1]
    return acf, overlap


def bright_orbital_features(image: np.ndarray) -> np.ndarray:
    """Represent elevated features by circular lobes around local maxima.

    Multiple maxima inside one connected bright region are retained as separate
    orbitals, so two touching lobes still count as two bright features.
    """
    smooth = ndimage.gaussian_filter(image, 0.7)
    bright_mask = smooth >= np.quantile(smooth, 0.88)
    maxima = (smooth == ndimage.maximum_filter(smooth, size=3)) & bright_mask
    labels, count = ndimage.label(bright_mask)
    centers = np.zeros_like(bright_mask, dtype=bool)
    for component in range(1, count + 1):
        region = labels == component
        component_maxima = maxima & region
        if np.any(component_maxima):
            centers |= component_maxima
        else:
            values = np.where(region, smooth, -np.inf)
            centers[np.unravel_index(np.argmax(values), values.shape)] = True
    yy, xx = np.mgrid[-2:3, -2:3]
    circular_orbital = (xx * xx + yy * yy) <= 4
    return ndimage.binary_dilation(centers, structure=circular_orbital)


def dark_orbital_features(image: np.ndarray) -> np.ndarray:
    """Represent depressed features with the same orbital construction."""
    return bright_orbital_features(-image)


def ordinary_orbital_features(image: np.ndarray) -> np.ndarray:
    """Represent morphology without preferring protrusions over depressions."""
    centered = image - float(np.median(image))
    return bright_orbital_features(np.abs(centered))


def has_distinct_light_structure(image: np.ndarray) -> bool:
    """Require a robust high-height excursion, not merely a top quantile."""
    median = float(np.median(image))
    robust_sigma = 1.4826 * float(np.median(np.abs(image - median))) + 1e-12
    return bool(float(np.quantile(image, 0.95)) >= median + 1.5 * robust_sigma)


def weighted_feature_f1(reference: np.ndarray, candidate: np.ndarray, f1_score) -> tuple[float, str]:
    """Apply 90/10 height weighting, with an 80% ordinary-feature fallback."""
    if has_distinct_light_structure(reference):
        score = (
            0.90 * f1_score(bright_orbital_features(reference), bright_orbital_features(candidate))
            + 0.10 * f1_score(dark_orbital_features(reference), dark_orbital_features(candidate))
        )
        return float(score), "90% high/light and 10% low/dark"
    score = (
        0.80 * f1_score(ordinary_orbital_features(reference), ordinary_orbital_features(candidate))
        + 0.20 * f1_score(dark_orbital_features(reference), dark_orbital_features(candidate))
    )
    return float(score), "80% ordinary morphology fallback and 20% low/dark"


def transformed_patch(image: np.ndarray, scale: float, rotation: float, mirrored: bool) -> np.ndarray:
    """Apply a small centered shape transform and restore the original array size."""
    transformed = image[:, ::-1] if mirrored else image
    transformed = ndimage.rotate(transformed, rotation, reshape=False, order=1, mode="nearest")
    if scale != 1.0:
        zoomed = ndimage.zoom(transformed, scale, order=1, mode="nearest")
        target_h, target_w = image.shape
        result = np.full((target_h, target_w), float(np.median(image)))
        src_y = max(0, (zoomed.shape[0] - target_h) // 2)
        src_x = max(0, (zoomed.shape[1] - target_w) // 2)
        dst_y = max(0, (target_h - zoomed.shape[0]) // 2)
        dst_x = max(0, (target_w - zoomed.shape[1]) // 2)
        height = min(target_h, zoomed.shape[0])
        width = min(target_w, zoomed.shape[1])
        result[dst_y:dst_y + height, dst_x:dst_x + width] = zoomed[src_y:src_y + height, src_x:src_x + width]
        transformed = result
    return transformed


def shift_metrics(image: np.ndarray, dx: int, dy: int) -> tuple[float, float, int, int, float, float, bool, float, float]:
    h, w = image.shape
    x0, x1 = max(0, -dx), min(w, w - dx)
    y0, y1 = max(0, -dy), min(h, h - dy)
    base = image[y0:y1, x0:x1]
    shifted = image[y0 + dy : y1 + dy, x0 + dx : x1 + dx]
    base_bright = bright_orbital_features(base)
    base_dark = dark_orbital_features(base)

    def f1_score(a: np.ndarray, b: np.ndarray) -> float:
        intersection = np.count_nonzero(a & b)
        precision = intersection / max(np.count_nonzero(a), 1)
        recall = intersection / max(np.count_nonzero(b), 1)
        return 2 * precision * recall / max(precision + recall, 1e-12)

    translation_corr = float(np.corrcoef(base.ravel(), shifted.ravel())[0, 1])
    translation_f1, _ = weighted_feature_f1(base, shifted, f1_score)

    correlation, f1 = translation_corr, translation_f1
    best_scale, best_rotation, best_mirrored = 1.0, 0.0, False
    # Limited tolerant comparison measures modest morphology changes without
    # allowing arbitrary transforms to manufacture unrelated matches.
    for scale in (0.90, 0.95, 1.0, 1.05, 1.10):
        for rotation in (-6.0, -3.0, 0.0, 3.0, 6.0):
            for mirrored in (False, True):
                if scale == 1.0 and rotation == 0.0 and not mirrored:
                    continue
                transformed = transformed_patch(shifted, scale, rotation, mirrored)
                corr = float(np.corrcoef(base.ravel(), transformed.ravel())[0, 1])
                transformed_f1, _ = weighted_feature_f1(base, transformed, f1_score)
                old_score = 0.55 * correlation + 0.45 * f1
                new_score = 0.55 * corr + 0.45 * transformed_f1
                complexity = 0.10 * abs(scale - 1.0) / 0.10 + 0.08 * abs(rotation) / 6.0 + 0.04 * mirrored
                if new_score - complexity > old_score:
                    correlation, f1 = corr, transformed_f1
                    best_scale, best_rotation, best_mirrored = scale, rotation, mirrored
    shifted_bright = bright_orbital_features(shifted)
    shifted_dark = dark_orbital_features(shifted)
    correlated_tiles = 0
    feature_tiles = 0
    for tile_y in range(4):
        for tile_x in range(4):
            ya, yb = tile_y * base.shape[0] // 4, (tile_y + 1) * base.shape[0] // 4
            xa, xb = tile_x * base.shape[1] // 4, (tile_x + 1) * base.shape[1] // 4
            local_corr = np.corrcoef(
                base[ya:yb, xa:xb].ravel(), shifted[ya:yb, xa:xb].ravel()
            )[0, 1]
            local_f1_bright = f1_score(
                base_bright[ya:yb, xa:xb], shifted_bright[ya:yb, xa:xb]
            )
            local_f1_dark = f1_score(
                base_dark[ya:yb, xa:xb], shifted_dark[ya:yb, xa:xb]
            )
            if has_distinct_light_structure(base[ya:yb, xa:xb]):
                local_f1 = 0.90 * local_f1_bright + 0.10 * local_f1_dark
            else:
                local_f1 = (
                    0.80 * f1_score(
                        ordinary_orbital_features(base[ya:yb, xa:xb]),
                        ordinary_orbital_features(shifted[ya:yb, xa:xb]),
                    )
                    + 0.20 * local_f1_dark
                )
            correlated_tiles += bool(local_corr > 0.20)
            feature_tiles += bool(local_f1 > 0.12)
    return (correlation, f1, correlated_tiles, feature_tiles, best_scale, best_rotation,
            best_mirrored, translation_corr, translation_f1)


def candidates(
    forward: np.ndarray, backward: np.ndarray, limit: int = 10
) -> tuple[list[Candidate], np.ndarray, np.ndarray]:
    af, overlap_f = normalized_acf(forward)
    ab, overlap_b = normalized_acf(backward)
    mean_acf = (af + ab) / 2
    overlap = np.minimum(overlap_f, overlap_b)
    h, w = forward.shape
    cy, cx = h - 1, w - 1
    # Use a compact neighborhood so a resolved secondary lobe close to the
    # origin is retained.  The origin itself is excluded below because it is
    # the zero-displacement self-correlation, not a second-tip separation.
    maxima = mean_acf == ndimage.maximum_filter(mean_acf, size=3)
    raw_peaks: list[tuple] = []

    for py, px in zip(*np.where(maxima)):
        dy, dx = py - cy, px - cx
        if dy < 0 or (dy == 0 and dx <= 0):
            continue  # retain one member of each automatically symmetric pair
        if math.hypot(dx, dy) < 4.5:
            continue  # exclude the zero-lag peak and its unresolved central core
        if abs(dx) > 0.40 * w or abs(dy) > 0.40 * h:
            continue
        overlap_fraction = overlap[py, px] / (h * w)
        if overlap_fraction < 0.50:
            continue

        local = mean_acf[max(0, py - 14) : py + 15, max(0, px - 14) : px + 15]
        median = np.median(local)
        robust_sigma = 1.4826 * np.median(np.abs(local - median)) + 1e-12
        prominence = (mean_acf[py, px] - median) / robust_sigma
        raw_peaks.append((float((af[py, px] + ab[py, px]) / 2), int(dx), int(dy),
                          float(af[py, px]), float(ab[py, px]), float(prominence),
                          float(overlap_fraction)))

    # A close doublet can appear only as a shoulder of the central ACF peak.
    # Add short vectors that are independently repeated by lobe geometry, but
    # never add a vector seen in only one or two bright components.
    morphology_image = (forward + backward) / 2
    for dx, dy, supporting_components in repeated_close_lobe_vectors(morphology_image):
        if math.hypot(dx, dy) < 4.5:
            continue
        if any(math.hypot(dx - item[1], dy - item[2]) <= 2.0 for item in raw_peaks):
            continue
        py, px = cy + dy, cx + dx
        if not (0 <= py < mean_acf.shape[0] and 0 <= px < mean_acf.shape[1]):
            continue
        overlap_fraction = float(overlap[py, px] / (h * w))
        if overlap_fraction < 0.50:
            continue
        local = mean_acf[max(0, py - 14):py + 15, max(0, px - 14):px + 15]
        median = float(np.median(local))
        robust_sigma = 1.4826 * float(np.median(np.abs(local - median))) + 1e-12
        prominence = float((mean_acf[py, px] - median) / robust_sigma)
        # The small bonus records independent morphology support without
        # allowing central-peak intensity alone to dominate candidate ranking.
        strength = float((af[py, px] + ab[py, px]) / 2 + 0.002 * supporting_components)
        raw_peaks.append((strength, int(dx), int(dy), float(af[py, px]),
                          float(ab[py, px]), prominence, overlap_fraction))

    raw_peaks.sort(reverse=True)
    found: list[Candidate] = []
    for _, dx, dy, acf_f, acf_b, prominence, overlap_fraction in raw_peaks[: max(limit * 3, limit)]:
        vector = np.asarray([dx, dy], dtype=float)
        norm = np.linalg.norm(vector)
        harmonic_count = 0
        if norm:
            unit = vector / norm
            for strength, other_dx, other_dy, *_ in raw_peaks:
                other = np.asarray([other_dx, other_dy], dtype=float)
                projection = float(other @ unit)
                perpendicular = float(np.linalg.norm(other - projection * unit))
                harmonic_error = abs(projection / norm - round(projection / norm))
                if strength > 0.10 and perpendicular <= 2.5 and harmonic_error <= 0.15:
                    harmonic_count += 1
        periodicity_penalty = min(1.0, max(0, harmonic_count - 1) / 3.0)
        mf = shift_metrics(forward, dx, dy)
        mb = shift_metrics(backward, dx, dy)
        found.append(
            Candidate(
                dx, dy, acf_f, acf_b, prominence,
                float(overlap_fraction), mf[0], mb[0], mf[1], mb[1],
                mf[2], mb[2], mf[3], mb[3], periodicity_penalty,
                mf[7], mb[7], mf[8], mb[8],
                mf[4], mb[4], mf[5], mb[5], mf[6], mb[6],
            )
        )
    found.sort(key=lambda item: (item.acf_forward + item.acf_backward) / 2, reverse=True)
    return found[:limit], af, ab


def classify(
    found: list[Candidate],
    detected_pairs: dict[str, list[dict]] | None = None,
    image_shape: tuple[int, int] | None = None,
) -> tuple[str, str, float]:
    detected_pairs = detected_pairs or {}

    def directional_family_count(items: list[Candidate], angle_tolerance_deg: float = 20.0) -> int:
        """Count sign-invariant vector orientations, not harmonics as extra tips."""
        family_angles: list[float] = []
        for item in items:
            angle = math.degrees(math.atan2(item.dy, item.dx)) % 180.0
            if not any(min(abs(angle - known), 180.0 - abs(angle - known)) <= angle_tolerance_deg
                       for known in family_angles):
                family_angles.append(angle)
        return len(family_angles)

    def pair_support(item: Candidate) -> tuple[int, int, float]:
        pairs = detected_pairs.get(f"{item.dx},{item.dy}", [])
        if not pairs:
            return 0, 0, 0.0
        regions: set[tuple[int, int]] = set()
        if image_shape:
            height, width = image_shape
            for pair in pairs:
                x, y = pair["feature_a_xy"]
                regions.add((min(3, 4 * y // max(height, 1)), min(3, 4 * x // max(width, 1))))
        return len(pairs), len(regions), float(np.mean([pair["confidence"] for pair in pairs]))

    def supported_directional_families(items: list[Candidate]) -> list[list[Candidate]]:
        """Aggregate nearby vectors before applying dense-field thresholds.

        A real separation can be split between adjacent pixel vectors. One
        member may carry the stronger ACF prominence while another carries
        more validated orbital pairs; treating them independently creates
        brittle threshold failures.
        """
        families: list[list[Candidate]] = []
        for item in items:
            angle = math.degrees(math.atan2(item.dy, item.dx)) % 180.0
            for family in families:
                known = math.degrees(math.atan2(family[0].dy, family[0].dx)) % 180.0
                if min(abs(angle - known), 180.0 - abs(angle - known)) <= 20.0:
                    family.append(item)
                    break
            else:
                families.append([item])

        supported: list[list[Candidate]] = []
        for family in families:
            supports = [pair_support(item) for item in family]
            total_pairs = sum(value[0] for value in supports)
            weighted_confidence = (
                sum(value[0] * value[2] for value in supports) / max(total_pairs, 1)
            )
            if (
                len(family) >= 2
                and total_pairs >= 40
                and max((item.acf_forward + item.acf_backward) / 2 for item in family) >= 0.15
                and max(item.prominence_z for item in family) >= 1.90
                and weighted_confidence >= 0.59
                and max(min(item.feature_tiles_forward, item.feature_tiles_backward) for item in family) >= 10
                and any(survives_lattice_check(item) for item in family)
            ):
                supported.append(family)
        return supported
    def effective_corr(item: Candidate) -> float:
        return max(
            min(item.translation_corr_forward, item.translation_corr_backward),
            min(item.corr_forward, item.corr_backward) - 0.05,
        )

    def effective_f1(item: Candidate) -> float:
        return max(
            min(item.translation_f1_forward, item.translation_f1_backward),
            min(item.f1_forward, item.f1_backward) - 0.04,
        )

    def independent_irregular_support(item: Candidate) -> bool:
        """Evidence beyond a normal repeating lattice at the same vector."""
        return (
            effective_corr(item) >= 0.30
            and effective_f1(item) >= 0.18
            and min(item.feature_tiles_forward, item.feature_tiles_backward) >= 5
        )

    def survives_lattice_check(item: Candidate) -> bool:
        return item.periodicity_penalty < 0.67 or independent_irregular_support(item)

    strong = [
        item for item in found
        if (
            item.prominence_z >= 3.0
            and effective_corr(item) >= 0.25
            and effective_f1(item) >= 0.20
            and min(item.feature_tiles_forward, item.feature_tiles_backward) >= 8
            and survives_lattice_check(item)
        )
        or (
            (item.acf_forward + item.acf_backward) / 2 >= 0.30
            and item.prominence_z >= 1.25
            and pair_support(item)[0] >= 5
            and pair_support(item)[1] >= 3
            and pair_support(item)[2] >= 0.55
            and survives_lattice_check(item)
        )
        or (
            # ACF prominence can be diluted by atomic-lattice peaks. Permit a
            # candidate only when direct irregular-feature evidence is strong
            # and broadly repeated; this is never an ACF-only decision.
            (item.acf_forward + item.acf_backward) / 2 >= 0.30
            and item.prominence_z >= 1.25
            and effective_corr(item) >= 0.35
            and effective_f1(item) >= 0.18
            and min(item.consistent_tiles_forward, item.consistent_tiles_backward) >= 12
            and min(item.feature_tiles_forward, item.feature_tiles_backward) >= 5
            and survives_lattice_check(item)
        )
        or (
            # Dense disordered fields can suppress whole-image correlation
            # even when the same touching orbital doublet is repeated many
            # times.  Accept a locally repeated family only when it has a
            # symmetric ACF counterpart, high prominence, many one-to-one
            # pairs, broad spatial coverage, and adequate mean confidence.
            (item.acf_forward + item.acf_backward) / 2 >= 0.10
            and item.prominence_z >= 2.0
            and pair_support(item)[0] >= 20
            and pair_support(item)[1] >= 6
            and pair_support(item)[2] >= 0.60
            and min(item.feature_tiles_forward, item.feature_tiles_backward) >= 6
            and survives_lattice_check(item)
        )
    ]
    moderate = [
        item for item in found
        if item.prominence_z >= 2.5
        and effective_corr(item) >= 0.18
            and effective_f1(item) >= 0.12
            and min(item.feature_tiles_forward, item.feature_tiles_backward) >= 5
            and survives_lattice_check(item)
    ]
    if found and all(item.periodicity_penalty >= 0.67 for item in found) and not any(
        independent_irregular_support(item) for item in found
    ):
        return "No clear evidence", (
            "All detected peak families are compatible with a normal lattice, and no independent "
            "elevated or irregular morphology repeats consistently at those displacements."
        ), 0.90
    if len(strong) >= 2 and directional_family_count(strong) >= 2:
        return "Likely multiple tip", "Several independently supported, distinctly oriented displacement families pass all evidence tests.", 0.90
    if len(strong) >= 2:
        return "Likely double tip", "Several supported vectors belong to one sign-invariant directional family and are treated as a changing or harmonic double-tip separation, not independent multiple tips.", 0.88
    if len(strong) == 1:
        return "Likely double tip", "One displacement vector passes all evidence tests.", 0.85
    family_supported = supported_directional_families(found)
    if len(family_supported) >= 2:
        return "Likely multiple tip", "Several aggregated directional families pass the dense-field orbital evidence tests.", 0.88
    if len(family_supported) == 1:
        return "Likely double tip", "Adjacent displacement vectors jointly support one repeated dense-field orbital family.", 0.86
    if moderate:
        return "Possible/uncertain", "At least one candidate passes moderate but not strong thresholds.", 0.55
    return "No clear evidence", "No autocorrelation candidate also reproduces irregular morphology consistently.", 0.80


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", type=Path, help="Nanonis .sxm or background-subtracted .npz file")
    parser.add_argument("--candidates", type=int, default=10, help="number of candidates to print")
    parser.add_argument("--json-output", type=Path, help="save downloadable evidence as JSON")
    parser.add_argument(
        "--acf-output", type=Path,
        help="save normalized forward, backward, and mean 2D autocorrelation maps as NPZ",
    )
    parser.add_argument(
        "--denoised-output", type=Path,
        help="save a background-corrected, noise-reduced grayscale STM image",
    )
    parser.add_argument("--flattened-output", type=Path, help="save the flattened/background-corrected STM image")
    parser.add_argument("--noise-mask-output", type=Path, help="save the isolated-glitch noise mask")
    parser.add_argument("--original-output", type=Path, help="save the original cropped STM image")
    parser.add_argument("--acf-image-output", type=Path, help="save an annotated global autocorrelation image")
    parser.add_argument("--feature-image-output", type=Path, help="save detected duplicated features marked on the STM image")
    args = parser.parse_args()

    z_forward, z_backward, pixel_x, pixel_y = load_input(args.input)
    original_shape = z_forward.shape
    z_forward, z_backward, crop = crop_common_finite(z_forward, z_backward)
    if args.original_output:
        args.original_output.parent.mkdir(parents=True, exist_ok=True)
        grayscale_image((z_forward + z_backward) / 2).save(args.original_output)
    flat_forward, reduced_forward, noise_forward, forward = preprocess_with_evidence(z_forward)
    flat_backward, reduced_backward, noise_backward, backward = preprocess_with_evidence(z_backward)
    if args.flattened_output:
        args.flattened_output.parent.mkdir(parents=True, exist_ok=True)
        grayscale_image((flat_forward + flat_backward) / 2).save(args.flattened_output)
    if args.denoised_output:
        save_noise_removed_image(reduced_forward, reduced_backward, args.denoised_output)
    if args.noise_mask_output:
        save_noise_mask(noise_forward | noise_backward, args.noise_mask_output)
    found, acf_forward, acf_backward = candidates(forward, backward, args.candidates)
    cleaned_mean = (forward + backward) / 2
    detected_pairs = {
        f"{item.dx},{item.dy}": feature_pairs(cleaned_mean, item.dx, item.dy)
        for item in found
    }
    label, reason, confidence = classify(found, detected_pairs, cleaned_mean.shape)
    if args.acf_image_output:
        save_acf_image((acf_forward + acf_backward) / 2, found, args.acf_image_output)
    if args.feature_image_output:
        def pair_map_rank(item: Candidate) -> tuple[float, float, int]:
            pairs = detected_pairs.get(f"{item.dx},{item.dy}", [])
            touching_support = sum(pair["confidence"] for pair in pairs if pair["touch_or_overlap"])
            mean_confidence = float(np.mean([pair["confidence"] for pair in pairs])) if pairs else 0.0
            return touching_support, mean_confidence, len(pairs)

        primary_item = max(found, key=pair_map_rank) if found else None
        primary_pairs = detected_pairs.get(f"{primary_item.dx},{primary_item.dy}", []) if primary_item else []
        save_feature_pair_image(cleaned_mean, primary_pairs, args.feature_image_output)

    if args.acf_output:
        args.acf_output.parent.mkdir(parents=True, exist_ok=True)
        h, w = forward.shape
        cy, cx = h - 1, w - 1
        np.savez_compressed(
            args.acf_output,
            acf_forward=acf_forward,
            acf_backward=acf_backward,
            acf_mean=(acf_forward + acf_backward) / 2,
            dx_pixels=np.arange(-(w - 1), w, dtype=int),
            dy_pixels=np.arange(-(h - 1), h, dtype=int),
            pixel_x_m=np.asarray(pixel_x),
            pixel_y_m=np.asarray(pixel_y),
            candidate_dx_pixels=np.asarray([item.dx for item in found], dtype=int),
            candidate_dy_pixels=np.asarray([item.dy for item in found], dtype=int),
            positive_peak_rows=np.asarray([cy + item.dy for item in found], dtype=int),
            positive_peak_columns=np.asarray([cx + item.dx for item in found], dtype=int),
            negative_peak_rows=np.asarray([cy - item.dy for item in found], dtype=int),
            negative_peak_columns=np.asarray([cx - item.dx for item in found], dtype=int),
        )

    if args.json_output:
        args.json_output.parent.mkdir(parents=True, exist_ok=True)
        payload = {
            "input": str(args.input),
            "classification": label,
            "confidence": confidence,
            "reason": reason,
            "coordinate_convention": "+x right, +y down; the negative autocorrelation counterpart is implicit",
            "pixel_size_nm": {"x": pixel_x * 1e9, "y": pixel_y * 1e9},
            "scan_shape": [int(z_forward.shape[0]), int(z_forward.shape[1])],
            "source_metadata": load_sxm_metadata(args.input),
            "forward_backward_morphology_correlation": float(np.corrcoef(forward.ravel(), backward.ravel())[0, 1]),
            "noise_filtering": {
                **NOISE_PARAMETERS,
                "forward_masked_pixels": int(np.count_nonzero(noise_forward)),
                "backward_masked_pixels": int(np.count_nonzero(noise_backward)),
                "union_masked_pixels": int(np.count_nonzero(noise_forward | noise_backward)),
                "union_masked_fraction": float(np.mean(noise_forward | noise_backward)),
            },
            "feature_weighting": "90% elevated bright orbital features and 10% depressed dark features; if no distinct light structure exists, 80% ordinary morphology and 20% dark features; detected peaks contribute equally regardless of height amplitude",
            "ordinary_analysis": {
                "forward_backward_correlation": float(np.corrcoef(forward.ravel(), backward.ravel())[0, 1]),
                "candidate_count": len(found),
            },
            "height_weighted_analysis": {
                "high_feature_weight": 1.0,
                "low_feature_weight": 0.1,
                "no_distinct_light_structure_weight": 0.8,
            },
            "feature_pair_map_vector": (
                {"dx": int(primary_item.dx), "dy": int(primary_item.dy)}
                if args.feature_image_output and primary_item is not None else None
            ),
            "candidates": [
                {
                    "dx": int(item.dx), "dy": int(item.dy),
                    "magnitude_pixels": float(math.hypot(item.dx, item.dy)),
                    "angle_degrees": float(math.degrees(math.atan2(item.dy, item.dx))),
                    **acf_peak_geometry((acf_forward + acf_backward) / 2, item.dx, item.dy),
                    "positive_peak_acf_location": {
                        "row": int((z_forward.shape[0] - 1) + item.dy),
                        "column": int((z_forward.shape[1] - 1) + item.dx),
                        "dx": int(item.dx), "dy": int(item.dy),
                    },
                    "negative_peak_acf_location": {
                        "row": int((z_forward.shape[0] - 1) - item.dy),
                        "column": int((z_forward.shape[1] - 1) - item.dx),
                        "dx": int(-item.dx), "dy": int(-item.dy),
                    },
                    "acf_forward": float(item.acf_forward), "acf_backward": float(item.acf_backward),
                    "prominence_z": float(item.prominence_z), "overlap_fraction": float(item.overlap_fraction),
                    "corr_forward": float(item.corr_forward), "corr_backward": float(item.corr_backward),
                    "f1_forward": float(item.f1_forward), "f1_backward": float(item.f1_backward),
                    "consistent_tiles_forward": int(item.consistent_tiles_forward),
                    "consistent_tiles_backward": int(item.consistent_tiles_backward),
                    "feature_tiles_forward": int(item.feature_tiles_forward),
                    "feature_tiles_backward": int(item.feature_tiles_backward),
                    "translation_only_morphology_forward": float(item.translation_corr_forward),
                    "translation_only_morphology_backward": float(item.translation_corr_backward),
                    "transformation_tolerant_morphology_forward": float(item.corr_forward),
                    "transformation_tolerant_morphology_backward": float(item.corr_backward),
                    "scale_ratio_forward": float(item.scale_forward),
                    "scale_ratio_backward": float(item.scale_backward),
                    "rotation_degrees_forward": float(item.rotation_forward_deg),
                    "rotation_degrees_backward": float(item.rotation_backward_deg),
                    "mirrored_forward": bool(item.mirrored_forward),
                    "mirrored_backward": bool(item.mirrored_backward),
                    "periodicity_penalty": float(item.periodicity_penalty),
                    "dx_nm": float(item.dx * pixel_x * 1e9),
                    "dy_nm": float(item.dy * pixel_y * 1e9),
                    "independent_regions_supporting": int(min(item.feature_tiles_forward, item.feature_tiles_backward)),
                    "spatial_consistency_score": float(min(item.feature_tiles_forward, item.feature_tiles_backward) / 16.0),
                    "lattice_correlation_present": bool(item.periodicity_penalty > 0),
                    "explained_only_by_normal_lattice": bool(item.periodicity_penalty >= 0.67 and min(item.feature_tiles_forward, item.feature_tiles_backward) < 5),
                    "non_lattice_duplicated_structures": bool(min(item.feature_tiles_forward, item.feature_tiles_backward) >= 5),
                    "duplicated_feature_pairs": detected_pairs.get(f"{item.dx},{item.dy}", []),
                }
                for item in found
            ],
        }
        args.json_output.write_text(json.dumps(payload, indent=2) + "\n")

    print(f"Classification: {label}")
    print(f"Confidence: {confidence:.2f}")
    print(f"Reason: {reason}")
    print(f"Scan: {z_forward.shape[1]} x {z_forward.shape[0]} pixels")
    if z_forward.shape != original_shape:
        print(f"Finite-data crop: x={crop[0]}:{crop[2]}, y={crop[1]}:{crop[3]} from {original_shape[1]} x {original_shape[0]}")
    print(f"Pixel calibration: x={pixel_x * 1e9:.6g} nm, y={pixel_y * 1e9:.6g} nm")
    print(f"Forward/backward morphology correlation: {np.corrcoef(forward.ravel(), backward.ravel())[0, 1]:.3f}")
    if args.json_output:
        print(f"Evidence JSON: {args.json_output}")
    if args.acf_output:
        print(f"Normalized 2D autocorrelation evidence: {args.acf_output}")
    if args.denoised_output:
        print(f"Noise-removed STM image: {args.denoised_output}")
    if args.flattened_output:
        print(f"Flattened STM image: {args.flattened_output}")
    if args.noise_mask_output:
        print(f"Noise mask: {args.noise_mask_output}")
    if args.original_output:
        print(f"Original STM image: {args.original_output}")
    if args.acf_image_output:
        print(f"Annotated autocorrelation image: {args.acf_image_output}")
    if args.feature_image_output:
        print(f"Marked duplicated-feature image: {args.feature_image_output}")
    print("\nCandidates (+x right, +y down; negative counterpart is implicit):")
    for item in found:
        cy, cx = z_forward.shape[0] - 1, z_forward.shape[1] - 1
        print(
            f"  ({item.dx:+d},{item.dy:+d}) px  "
            f"({item.dx * pixel_x * 1e9:+.3f},{item.dy * pixel_y * 1e9:+.3f}) nm  "
            f"ACF F/B={item.acf_forward:.3f}/{item.acf_backward:.3f}  "
            f"prominence={item.prominence_z:.2f} sigma  "
            f"corr F/B={item.corr_forward:.3f}/{item.corr_backward:.3f}  "
            f"feature-F1 F/B={item.f1_forward:.3f}/{item.f1_backward:.3f}  "
            f"correlated tiles F/B={item.consistent_tiles_forward}/16/{item.consistent_tiles_backward}/16  "
            f"feature tiles F/B={item.feature_tiles_forward}/16/{item.feature_tiles_backward}/16"
            f"  periodicity penalty={item.periodicity_penalty:.2f}"
            f"  ACF pair locations=(r{cy + item.dy},c{cx + item.dx})/"
            f"(r{cy - item.dy},c{cx - item.dx})"
        )


if __name__ == "__main__":
    main()
