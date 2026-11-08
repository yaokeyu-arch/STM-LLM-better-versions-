---
name: stm-double-tip-analysis
description: Analyze Nanonis STM/SXM images for evidence of double or multiple tips using conservative preprocessing, normalized autocorrelation, one-to-one morphology pairing, sliced consistency, and lattice rejection. Use for individual scans, folders of scans, reruns, and evidence or evaluation exports.
---

# STM Double-Tip Analysis

Use the scripts packaged with this skill. Apply `$saving-protocol` for every run so overall results and per-image evidence remain organized and earlier runs are not overwritten.

## Scripts

- `scripts/main.py`: preferred skill-only entry point for one SXM file or a folder; applies the saving protocol, archives source data and exact materials, runs the packaged analysis scripts, and writes a checksummed inventory.
- `scripts/stm_background_subtract.py`: conservative plane and scan-line correction.
- `scripts/stm_double_tip_analysis.py`: whole-image autocorrelation and morphology analysis.
- `scripts/Sliced Autocorrection/sliced_autocorrelation.py`: acquisition-ordered sliding-window analysis.
- `scripts/run_stm_double_tip_pipeline.py`: complete single-image evidence pipeline.
- `scripts/batch_analyze_folder.py`: complete folder analysis using the saving protocol.
- `scripts/build_stm_evaluation_workbooks.mjs`: ground-truth comparison workbooks.
- `scripts/build_run_summary_workbook.mjs`: overall interpretation workbook.

## Required interpretation principles

- Preserve original height data and save original, flattened, and conservatively noise-reduced versions.
- Never classify from an autocorrelation peak alone. Require real morphology duplication and spatial or sliced support.
- Treat touching orbital lobes as distinct peaks when local geometry supports separation.
- Detect lobes at multiple conservative scales: use a lightly smoothed map to retain weaker shoulders and an approximately top-20% bright support region, together with a robust intensity floor. Do not require both members of a visible orbital pair to fall within the top 12% of pixels.
- A weaker lobe beside a stronger maximum remains a candidate when it is locally maximal, lies in or immediately beside the same bright support region, and subsequently passes the morphology, height, displacement, and one-to-one checks.
- Do not require two separate discrete local maxima for a touching pair. If a compact bright core is strongly elongated, deblend it along its principal axis and retain two candidate lobe centres even when the saddle is shallow. The candidates must still pass the downstream displacement, morphology, one-to-one, and spatial-consistency checks.
- Pair individual lobes one-to-one. Never reuse a peak, join two composite doublets as one pair, or pair component centroids when individual lobes are detectable.
- Resolve partner assignments globally across all plausible displacement families. Do not finalize pairs independently for each ACF vector: a locally good short-vector match must not take a lobe away from a better long-vector orbital match.
- When competing partners exist, compare the normalized shape and contour of the bright orbital core, not mainly the surrounding patch or absolute peak height. Two lobes may be a valid pair even when their brightness differs substantially.
- Prefer the assignment that preserves the complete repeated orbital arrangement and yields the same displacement for other independent structures. Nearest-neighbor distance is not, by itself, a pairing criterion.
- Preserve user-verified pair corrections in the run evidence and use them as regression cases. For `12_23_21_017`, the verified pair `(188,53) -> (238,55)` has displacement `(+50,+2) px`; the former short-vector partners for those lobes are incorrect.
- Search for short, repeated intra-component lobe vectors even when they appear only as a shoulder of the central autocorrelation peak. Require the short vector in at least three independent bright components before adding it as a candidate.
- Exclude the unresolved central ACF core: candidate magnitude must be at least `4.5 px`. A short vector outside that core may be retained only when repeated lobe geometry independently supports it.
- Apply the same `4.5 px` unresolved-core exclusion independently in every sliced/windowed autocorrelation. A `(0,1)`, `(1,1)`, or similarly tiny shift is central self-correlation and must never create partial double-tip evidence.
- Validate every displayed feature pair before it contributes evidence. Require translation similarity of at least `0.30`, transformation-tolerant similarity of at least `0.35`, height similarity of at least `0.45`, and shape similarity of at least `0.50`, unless both lobes belong to the same resolved touching component.
- Scale displacement tolerance with vector length. For short vectors, use approximately `max(1.25 px, 20% of vector magnitude)`, capped at `3 px`; do not use a fixed ±3-pixel tolerance for a five-pixel lobe separation.
- Reject inter-doublet connections when the proposed members have more than `0.50` deformation, strongly unequal component sizes, or poor height similarity. Do not retain a visually convenient line merely because its displacement resembles an ACF vector.
- Number pair labels only after invalid pairs have been removed, so every `Pxx` label in the image corresponds exactly to the same accepted record in the evidence JSON.
- Display every accepted pair for the selected vector up to the documented output limit; do not display only a hand-picked example. Record the total accepted-pair count so omitted display labels are explicit if the limit is reached.
- Count validated feature pairs before making the final whole-image classification. At least five accepted pairs spanning at least three independent 4×4 spatial regions, with mean pair confidence of at least `0.55`, can provide strong morphology evidence when the associated ACF pair is meaningful and the lattice-only test is passed.
- For a dense disordered field, do not require high whole-image shifted correlation when unrelated surrounding morphology dilutes it. A vector may still provide strong evidence when it has a symmetric ACF pair of at least about `0.10`, prominence of at least `2 sigma`, at least 20 validated one-to-one orbital pairs across at least six independent regions, mean pair confidence of at least `0.60`, and broad feature-tile support. This exception must not be applied to a normal lattice-only family.
- Aggregate neighboring, sign-invariant displacement vectors before applying dense-field thresholds. Do not reject a family merely because one adjacent pixel vector contains the stronger ACF prominence while another contains the larger and more confident set of validated pairs. For an aggregate exception, require at least two vectors within about `20°`, at least 40 total validated pairs, ACF strength of at least `0.15`, maximum prominence around `1.9 sigma`, pair-weighted confidence of at least `0.59`, broad feature-tile support, and passage of the lattice-only rejection.
- Select the displayed feature-pair map by validated touching-lobe support and pair confidence, not merely by the highest raw ACF value. Record the vector used for the feature map in the evidence JSON.
- If several short lattice-like vector families compete and no single non-lattice vector is clearly supported, return `Possible/uncertain` and do not present any one family as confirmed ground truth.
- Give high/topographic features substantially more weight while retaining low-feature evidence at 10% weight.
- Reject normal lattice periodicity unless lattice sites themselves or non-lattice features show consistent duplicated morphology.
- Allow only limited, reported scale, rotation, deformation, and mirror transformations.
- Accept partial-image evidence and smoothly evolving vectors when morphology remains consistent.
- Distinguish `Likely double tip`, `Likely multiple tip`, `Possible/uncertain`, and `No clear evidence`.
- Do not call several collinear, sign-reversed, harmonic, or smoothly varying magnitudes separate tips. Cluster supported vectors using sign-invariant orientation; require at least two independently supported directional families separated by more than approximately `20°` before returning `Likely multiple tip`. Otherwise report a changing or harmonic double-tip vector family.

## Running a group

Run `scripts/batch_analyze_folder.py <image-group-folder>`. By default it creates a non-overwriting run under the workspace `runs/` folder, named after the image-group folder. Each image receives its own evidence folder and group-level results go in `overall-results/`.

Do not read ground-truth labels until blind predictions have been written. Use ground truth only for subsequent evaluation and error analysis.
