#!/usr/bin/env python3
"""Detect persistent single-tip -> double-tip transitions in Nanonis SXM scans.

The script follows the scan's metadata-defined slow-scan acquisition order,
uses whole-image autocorrelation only to propose reference displacements, and
then evaluates overlapping windows using autocorrelation plus direct matching
of irregular morphology. A transition is reported only when double-tip
evidence persists through every remaining window.

Example:
    python sliced_autocorrelation.py scan.sxm --output-dir results

Dependencies: NumPy and SciPy.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import re
from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np
from PIL import Image, ImageDraw
from scipy import ndimage


@dataclass
class Peak:
    dx: int
    dy: int
    strength: float
    prominence: float
    symmetry: float


@dataclass
class WindowResult:
    index: int
    acquisition_start: int
    acquisition_end: int
    image_row_start: int
    image_row_end: int
    dx_px: int | None
    dy_px: int | None
    dx_nm: float | None
    dy_nm: float | None
    peak_strength: float
    peak_prominence: float
    pair_symmetry: float
    morphology_correlation: float
    irregular_feature_f1: float
    vector_consistency: float
    neighbor_agreement: float
    periodicity_penalty: float
    evidence_score: float
    scale_factor: float = 1.0
    size_change_percent: float = 0.0
    rotation_degrees: float = 0.0
    mirrored: bool = False
    translation_only_correlation: float = 0.0
    translation_only_feature_f1: float = 0.0


def header_field(header: str, name: str) -> str:
    match = re.search(rf":{re.escape(name)}:\s*\n([^\n]+)", header)
    if not match:
        raise ValueError(f"Missing SXM header field: {name}")
    return match.group(1).strip()


def load_sxm(path: Path) -> tuple[np.ndarray, np.ndarray, dict]:
    raw = path.read_bytes()
    marker = b":SCANIT_END:"
    marker_at = raw.find(marker)
    if marker_at < 0:
        raise ValueError("Not a recognized Nanonis SXM file")
    header = raw[: marker_at + len(marker)].decode("latin1", "replace")
    nx, ny = map(int, header_field(header, "SCAN_PIXELS").split()[:2])
    range_x, range_y = map(float, header_field(header, "SCAN_RANGE").split()[:2])
    scan_direction = header_field(header, "SCAN_DIR").strip().lower()

    table = re.search(r":DATA_INFO:\s*\n(.*?)(?=\n:SCANIT_END:)", header, re.S)
    if not table:
        raise ValueError("Missing DATA_INFO table")
    rows = [line.strip() for line in table.group(1).splitlines() if line.strip()]
    channels: list[tuple[str, str]] = []
    for row in rows[1:]:
        columns = row.split("\t")
        if len(columns) >= 4:
            channels.append((columns[1].strip(), columns[3].strip().lower()))
    z_index = next((i for i, item in enumerate(channels) if item[0] == "Z"), None)
    if z_index is None:
        raise ValueError("No Z channel found")

    frame_counts = [2 if direction == "both" else 1 for _, direction in channels]
    frame_count = sum(frame_counts)
    data_offset = len(raw) - frame_count * nx * ny * 4
    frames = np.frombuffer(raw, dtype=">f4", offset=data_offset).reshape(frame_count, ny, nx)
    first = sum(frame_counts[:z_index])
    forward = frames[first].astype(float)
    backward = (
        frames[first + 1][:, ::-1].astype(float)
        if frame_counts[z_index] == 2
        else forward.copy()
    )
    metadata = {
        "scan_direction": scan_direction,
        "pixel_x_m": range_x / nx,
        "pixel_y_m": range_y / ny,
        "original_width": nx,
        "original_height": ny,
    }
    return forward, backward, metadata


def longest_true_run(mask: np.ndarray) -> tuple[int, int]:
    best_start = best_end = start = 0
    active = False
    for index, value in enumerate(np.r_[mask, False]):
        if value and not active:
            start, active = index, True
        elif not value and active:
            if index - start > best_end - best_start:
                best_start, best_end = start, index
            active = False
    return best_start, best_end


def finite_crop(
    forward: np.ndarray, backward: np.ndarray
) -> tuple[np.ndarray, np.ndarray, tuple[int, int, int, int]]:
    valid = np.isfinite(forward) & np.isfinite(backward)
    y0, y1 = longest_true_run(np.all(valid, axis=1))
    if y1 <= y0:
        raise ValueError("No common complete scan lines")
    x0, x1 = longest_true_run(np.all(valid[y0:y1], axis=0))
    if x1 <= x0 or (x1 - x0) * (y1 - y0) < 0.20 * valid.size:
        raise ValueError("Too little common finite data for reliable analysis")
    return forward[y0:y1, x0:x1], backward[y0:y1, x0:x1], (x0, y0, x1, y1)


def subtract_background(z: np.ndarray) -> np.ndarray:
    """Conservative plane and isolated-line-offset correction."""
    h, w = z.shape
    y, x = np.mgrid[:h, :w]
    design = np.c_[np.ones(z.size), x.ravel(), y.ravel()]
    plane = design @ np.linalg.lstsq(design, z.ravel(), rcond=None)[0]
    residual = z - plane.reshape(h, w)
    row_medians = np.median(residual, axis=1)
    filter_size = max(5, h // 15 | 1)
    line_offsets = row_medians - ndimage.median_filter(row_medians, filter_size)
    return residual - line_offsets[:, None]


def morphology_image(z: np.ndarray) -> np.ndarray:
    """Minimal denoising; broad subtraction is only for candidate comparison."""
    h, w = z.shape
    image = ndimage.gaussian_filter(z, 0.65)
    image -= ndimage.gaussian_filter(z, max(6.0, min(h, w) / 5.0))
    scale = image.std()
    return (image - image.mean()) / scale if scale else image * 0


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
    center = acf[h - 1, w - 1]
    if center:
        acf /= center
    return acf, overlap


def find_peaks(
    image: np.ndarray,
    max_peaks: int = 12,
    min_radius: float = 4.5,
    max_fraction: float = 0.40,
) -> list[Peak]:
    acf, overlap = normalized_acf(image)
    h, w = image.shape
    cy, cx = h - 1, w - 1
    # Exclude the unresolved zero-lag core. Correlation at one or two pixels is
    # dominated by the central peak and cannot establish a second-tip vector.
    maxima = acf == ndimage.maximum_filter(acf, size=3)
    peaks: list[Peak] = []
    for py, px in zip(*np.where(maxima)):
        dy, dx = py - cy, px - cx
        if dy < 0 or (dy == 0 and dx <= 0):
            continue
        if np.hypot(dx, dy) < min_radius:
            continue
        if abs(dx) > max_fraction * w or abs(dy) > max_fraction * h:
            continue
        if overlap[py, px] < 0.50 * h * w:
            continue
        local = acf[max(0, py - 12) : py + 13, max(0, px - 12) : px + 13]
        median = float(np.median(local))
        robust_sigma = 1.4826 * float(np.median(np.abs(local - median))) + 1e-12
        positive = float(acf[py, px])
        negative = float(acf[cy - dy, cx - dx])
        symmetry = min(abs(positive), abs(negative)) / max(abs(positive), abs(negative), 1e-12)
        peaks.append(Peak(int(dx), int(dy), float(positive),
                          float((positive - median) / robust_sigma), float(symmetry)))
    peaks.sort(key=lambda peak: peak.strength, reverse=True)
    return peaks[:max_peaks]


def f1_score(a: np.ndarray, b: np.ndarray) -> float:
    intersection = np.count_nonzero(a & b)
    precision = intersection / max(np.count_nonzero(a), 1)
    recall = intersection / max(np.count_nonzero(b), 1)
    return 2 * precision * recall / max(precision + recall, 1e-12)


def bright_orbital_features(image: np.ndarray) -> np.ndarray:
    """Create circular bright-feature groups while retaining touching lobes."""
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
    """Create matching orbital groups for depressed/dark features."""
    return bright_orbital_features(-image)


def ordinary_orbital_features(image: np.ndarray) -> np.ndarray:
    centered = image - float(np.median(image))
    return bright_orbital_features(np.abs(centered))


def has_distinct_light_structure(image: np.ndarray) -> bool:
    median = float(np.median(image))
    robust_sigma = 1.4826 * float(np.median(np.abs(image - median))) + 1e-12
    return bool(float(np.quantile(image, 0.95)) >= median + 1.5 * robust_sigma)


def weighted_feature_f1(reference: np.ndarray, candidate: np.ndarray) -> tuple[float, str]:
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


def morphology_match(image: np.ndarray, dx: int, dy: int) -> tuple[float, float, float, float, bool, float, float]:
    h, w = image.shape
    x0, x1 = max(0, -dx), min(w, w - dx)
    y0, y1 = max(0, -dy), min(h, h - dy)
    if x1 - x0 < 8 or y1 - y0 < 8:
        return 0.0, 0.0, 1.0, 0.0, False, 0.0, 0.0
    base = image[y0:y1, x0:x1]
    shifted = image[y0 + dy : y1 + dy, x0 + dx : x1 + dx]
    bright_a = bright_orbital_features(base)
    dark_a = dark_orbital_features(base)
    translation_correlation = float(np.corrcoef(base.ravel(), shifted.ravel())[0, 1])
    translation_feature_f1, _ = weighted_feature_f1(base, shifted)
    best_correlation = translation_correlation
    best_feature_f1 = translation_feature_f1
    best_scale, best_rotation, best_mirrored = 1.0, 0.0, False
    for scale in (0.90, 0.95, 1.0, 1.05, 1.10):
        for rotation in (-6.0, -3.0, 0.0, 3.0, 6.0):
            for mirrored in (False, True):
                if scale == 1.0 and rotation == 0.0 and not mirrored:
                    continue
                transformed = transformed_patch(shifted, scale, rotation, mirrored)
                correlation = float(np.corrcoef(base.ravel(), transformed.ravel())[0, 1])
                feature_f1, _ = weighted_feature_f1(base, transformed)
                old_score = 0.55 * best_correlation + 0.45 * best_feature_f1
                new_score = 0.55 * correlation + 0.45 * feature_f1
                complexity = 0.10 * abs(scale - 1.0) / 0.10 + 0.08 * abs(rotation) / 6.0 + 0.04 * mirrored
                if new_score - complexity > old_score:
                    best_correlation, best_feature_f1 = correlation, feature_f1
                    best_scale, best_rotation, best_mirrored = scale, rotation, mirrored
    return (best_correlation, best_feature_f1, best_scale, best_rotation, best_mirrored,
            translation_correlation, translation_feature_f1)


def vector_similarity(dx: int, dy: int, references: list[Peak], tolerance: float = 4.0) -> float:
    if not references:
        return 0.0
    distance = min(math.hypot(dx - peak.dx, dy - peak.dy) for peak in references)
    return float(math.exp(-0.5 * (distance / tolerance) ** 2))


def periodicity_penalty(peaks: list[Peak], chosen: Peak) -> float:
    """Penalize comb-like peak families compatible with regular periodicity."""
    vector = np.array([chosen.dx, chosen.dy], dtype=float)
    norm = np.linalg.norm(vector)
    if norm == 0:
        return 1.0
    unit = vector / norm
    harmonic_count = 0
    for peak in peaks:
        candidate = np.array([peak.dx, peak.dy], dtype=float)
        projection = float(candidate @ unit)
        perpendicular = float(np.linalg.norm(candidate - projection * unit))
        harmonic = abs(projection / norm - round(projection / norm))
        if perpendicular <= 2.5 and harmonic <= 0.15 and peak.strength > 0.10:
            harmonic_count += 1
    return min(1.0, max(0, harmonic_count - 1) / 3.0)


def clamp01(value: float) -> float:
    return float(np.clip(value, 0.0, 1.0))


def score_windows(
    image: np.ndarray,
    windows: list[tuple[int, int, int, int]],
    global_peaks: list[Peak],
    pixel_x_m: float,
    pixel_y_m: float,
) -> list[WindowResult]:
    provisional: list[dict] = []
    for index, (acq_start, acq_end, row_start, row_end) in enumerate(windows):
        window = image[row_start:row_end]
        peaks = find_peaks(window)
        best = None
        for peak in peaks:
            corr, feature_f1, scale, rotation, mirrored, translation_corr, translation_f1 = morphology_match(window, peak.dx, peak.dy)
            scoring_corr = max(translation_corr, corr - 0.05)
            scoring_f1 = max(translation_f1, feature_f1 - 0.04)
            global_consistency = vector_similarity(peak.dx, peak.dy, global_peaks)
            penalty = periodicity_penalty(peaks, peak)
            # Symmetry is recorded but lightly weighted: real-image ACF symmetry
            # is mathematical and therefore cannot independently prove double tip.
            # A lattice-family penalty is reduced when elevated irregular
            # morphology independently aligns at the same displacement.
            independent_irregular_support = scoring_corr >= 0.30 and scoring_f1 >= 0.18
            effective_penalty = penalty * (0.25 if independent_irregular_support else 1.0)
            score = (
                0.16 * clamp01(peak.strength / 0.45)
                + 0.10 * clamp01(peak.prominence / 3.0)
                + 0.04 * peak.symmetry
                + 0.25 * clamp01(scoring_corr / 0.50)
                + 0.30 * clamp01(scoring_f1 / 0.20)
                + 0.15 * global_consistency
                - 0.18 * effective_penalty
            )
            candidate = (score, peak, corr, feature_f1, global_consistency, effective_penalty,
                         scale, rotation, mirrored)
            candidate += (translation_corr, translation_f1)
            if best is None or candidate[0] > best[0]:
                best = candidate
        provisional.append({"best": best, "index": index, "bounds": (acq_start, acq_end, row_start, row_end)})

    results: list[WindowResult] = []
    for item in provisional:
        best = item["best"]
        acq_start, acq_end, row_start, row_end = item["bounds"]
        if best is None:
            results.append(WindowResult(item["index"], acq_start, acq_end, row_start, row_end,
                                        None, None, None, None, 0, 0, 0, 0, 0, 0, 0, 1, 0))
            continue
        base_score, peak, corr, feature_f1, consistency, penalty, scale, rotation, mirrored, translation_corr, translation_f1 = best
        neighbors: list[float] = []
        for neighbor_index in (item["index"] - 1, item["index"] + 1):
            if 0 <= neighbor_index < len(provisional) and provisional[neighbor_index]["best"]:
                neighbor_peak = provisional[neighbor_index]["best"][1]
                distance = math.hypot(peak.dx - neighbor_peak.dx, peak.dy - neighbor_peak.dy)
                neighbors.append(math.exp(-0.5 * (distance / 4.0) ** 2))
        agreement = float(np.mean(neighbors)) if neighbors else 0.0
        score = clamp01(0.85 * base_score + 0.15 * agreement)
        results.append(WindowResult(
            item["index"], acq_start, acq_end, row_start, row_end,
            peak.dx, peak.dy, peak.dx * pixel_x_m * 1e9, peak.dy * pixel_y_m * 1e9,
            peak.strength, peak.prominence, peak.symmetry, corr, feature_f1,
            consistency, agreement, penalty, score,
            scale, (scale - 1.0) * 100, rotation, mirrored,
            translation_corr, translation_f1,
        ))
    return results


def make_windows(height: int, fraction: float, overlap: float, scan_direction: str) -> list[tuple[int, int, int, int]]:
    size = max(16, min(height, round(height * fraction)))
    step = max(1, round(size * (1.0 - overlap)))
    starts = list(range(0, max(1, height - size + 1), step))
    if not starts or starts[-1] != height - size:
        starts.append(height - size)
    row_windows = [(start, start + size) for start in starts]
    # Image rows are stored top-to-bottom; SCAN_DIR defines temporal order.
    if scan_direction.startswith("up"):
        row_windows.reverse()
    return [(i * step, i * step + size, start, end) for i, (start, end) in enumerate(row_windows)]


def persistent_transition(results: list[WindowResult], threshold: float, minimum_tail: int) -> tuple[int | None, str]:
    scores = np.asarray([result.evidence_score for result in results])
    valid_changes: list[tuple[float, int]] = []
    for change in range(1, len(scores) - minimum_tail + 1):
        before, after = scores[:change], scores[change:]
        # Fundamental rule: every remaining window must retain double-tip evidence.
        if np.all(after >= threshold) and np.all(before < threshold):
            separation = float(np.mean(after) - np.mean(before))
            valid_changes.append((separation, change))
    if not valid_changes:
        if len(scores) >= minimum_tail and np.all(scores >= threshold):
            return 0, "Double-tip evidence is present from the beginning and persists to the end."
        return None, "No single change point produces persistent evidence through the end."
    _, change = max(valid_changes)
    return change, "Evidence crosses the threshold and remains above it through the final window."


def eligible_correlation_windows(results: list[WindowResult], threshold: float) -> list[WindowResult]:
    """Return windows with direct morphology-supported double-tip evidence."""
    return [
        result for result in results
        if result.dx_px is not None
        and result.evidence_score >= threshold
        and max(result.translation_only_correlation, result.morphology_correlation - 0.05) >= 0.18
        and max(result.translation_only_feature_f1, result.irregular_feature_f1 - 0.04) >= 0.12
    ]


def partial_correlation_windows(results: list[WindowResult], threshold: float) -> list[int]:
    """Return windows participating in a consistent morphology-supported partial match.

    A partial result requires at least two windows with direct morphology and
    bright-feature overlap at approximately the same displacement.  This keeps
    a lone high autocorrelation window from becoming a double-tip decision.
    """
    eligible = eligible_correlation_windows(results, threshold)
    matched: set[int] = set()
    for i, first in enumerate(eligible):
        for second in eligible[i + 1:]:
            distance = math.hypot(first.dx_px - second.dx_px, first.dy_px - second.dy_px)
            if distance <= 4.0:
                matched.update((first.index, second.index))
    return sorted(matched)


def distributed_changing_vector_windows(results: list[WindowResult], threshold: float) -> list[int]:
    """Detect broad sliced evidence whose displacement changes across the scan."""
    eligible = eligible_correlation_windows(results, threshold)
    if len(eligible) < max(3, math.ceil(0.50 * len(results))):
        return []
    indices = sorted(result.index for result in eligible)
    last = max(len(results) - 1, 1)
    if indices[0] > last / 3 or indices[-1] < 2 * last / 3:
        return []
    vectors = [(result.dx_px, result.dy_px) for result in eligible]
    maximum_change = max(
        math.hypot(ax - bx, ay - by)
        for i, (ax, ay) in enumerate(vectors)
        for bx, by in vectors[i + 1:]
    )
    return indices if maximum_change > 4.0 else []


def analyze(args: argparse.Namespace) -> tuple[dict, list[WindowResult]]:
    forward, backward, metadata = load_sxm(args.sxm)
    forward, backward, crop = finite_crop(forward, backward)
    forward = subtract_background(forward)
    backward = subtract_background(backward)
    # Average independent traces after registration to reduce line-direction noise.
    combined = morphology_image((forward + backward) / 2)
    global_peaks = find_peaks(combined)
    windows = make_windows(combined.shape[0], args.window_fraction, args.overlap, metadata["scan_direction"])
    results = score_windows(
        combined, windows, global_peaks, metadata["pixel_x_m"], metadata["pixel_y_m"]
    )
    change, explanation = persistent_transition(results, args.threshold, args.minimum_tail)
    partial_windows = partial_correlation_windows(results, args.threshold)
    changing_windows = distributed_changing_vector_windows(results, args.threshold)
    if change is None:
        if changing_windows:
            classification = "Partial double-tip correlations with changing displacement"
            explanation = (
                "Duplicated morphology is distributed across the scan while the displacement "
                f"changes between windows {', '.join(map(str, changing_windows))}."
            )
        elif partial_windows:
            classification = "Partial double-tip correlations"
            explanation = (
                "Morphology-supported correlations at a consistent displacement occur in "
                f"windows {', '.join(map(str, partial_windows))}, but do not persist to scan end."
            )
        else:
            classification = "No persistent transition"
    elif change == 0:
        classification = "Persistent double-tip evidence from scan start"
    else:
        classification = "Persistent single-to-double-tip transition"
    summary = {
        "input": str(args.sxm),
        "classification": classification,
        "explanation": explanation,
        "scan_direction_metadata": metadata["scan_direction"],
        "acquisition_window_order": "bottom-to-top" if metadata["scan_direction"].startswith("up") else "top-to-bottom",
        "finite_crop_xyxy": crop,
        "analyzed_shape": list(combined.shape),
        "pixel_size_nm": [metadata["pixel_x_m"] * 1e9, metadata["pixel_y_m"] * 1e9],
        "window_fraction": args.window_fraction,
        "overlap": args.overlap,
        "evidence_threshold": args.threshold,
        "feature_weighting": "90% elevated bright orbital features and 10% depressed dark features; if no distinct light structure exists, 80% ordinary morphology and 20% dark features; detected peaks contribute equally regardless of height amplitude",
        "change_window_index": change,
        "partial_correlation_window_indices": partial_windows,
        "changing_displacement_window_indices": changing_windows,
        "vector_evolution": [
            {"window": result.index, "dx": result.dx_px, "dy": result.dy_px,
             "score": result.evidence_score}
            for result in results if result.dx_px is not None
        ],
        "global_reference_peaks": [asdict(peak) for peak in global_peaks],
        "scores": [result.evidence_score for result in results],
    }
    return summary, results


def save_results(output_dir: Path, stem: str, summary: dict, results: list[WindowResult]) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    json_path = output_dir / f"{stem}_sliced_summary.json"
    csv_path = output_dir / f"{stem}_window_scores.csv"
    json_path.write_text(json.dumps(summary, indent=2) + "\n")
    with csv_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(asdict(results[0]).keys()))
        writer.writeheader()
        writer.writerows(asdict(result) for result in results)
    height, width = summary["analyzed_shape"]
    scale = min(1.0, 900 / max(width, height))
    canvas = Image.new("RGB", (max(1, round(width * scale)), max(1, round(height * scale))), (30, 30, 30))
    draw = ImageDraw.Draw(canvas, "RGBA")
    for result in results:
        y0, y1 = round(result.image_row_start * scale), round(result.image_row_end * scale)
        intensity = int(255 * float(np.clip(result.evidence_score, 0, 1)))
        color = (255, 60, 30, 65) if result.evidence_score >= summary["evidence_threshold"] else (30, 140, 255, 45)
        draw.rectangle((0, y0, canvas.width - 1, max(y0, y1 - 1)), outline=color[:3] + (210,), fill=color)
        label = f"W{result.index} d=({result.dx_px},{result.dy_px}) S={result.evidence_score:.2f}"
        draw.text((4, y0 + 2), label, fill=(255, 255, intensity, 255))
    window_image = output_dir / f"{stem}_sliced_window_map.png"
    canvas.save(window_image)
    print(f"Classification: {summary['classification']}")
    print(f"Acquisition order: {summary['acquisition_window_order']} ({summary['scan_direction_metadata']})")
    print(f"Explanation: {summary['explanation']}")
    print(f"Scores: {', '.join(f'{score:.3f}' for score in summary['scores'])}")
    print(f"Summary: {json_path}")
    print(f"Window measurements: {csv_path}")
    print(f"Sliced evidence map: {window_image}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("sxm", type=Path)
    parser.add_argument("--window-fraction", type=float, default=0.20,
                        help="window height as fraction of valid image height (default: 0.20)")
    parser.add_argument("--overlap", type=float, default=0.625,
                        help="fractional overlap between windows (default: 0.625)")
    parser.add_argument("--threshold", type=float, default=0.62,
                        help="per-window double-tip evidence threshold (default: 0.62)")
    parser.add_argument("--minimum-tail", type=int, default=3,
                        help="minimum persistent post-transition windows (default: 3)")
    parser.add_argument("--output-dir", type=Path, default=Path.cwd())
    args = parser.parse_args()
    if not 0.15 <= args.window_fraction <= 0.25:
        parser.error("--window-fraction must be between 0.15 and 0.25")
    if not 0.50 <= args.overlap <= 0.75:
        parser.error("--overlap must be between 0.50 and 0.75")
    if not 0 < args.threshold < 1:
        parser.error("--threshold must be between 0 and 1")
    summary, results = analyze(args)
    save_results(args.output_dir, args.sxm.stem, summary, results)


if __name__ == "__main__":
    main()
