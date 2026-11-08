---
name: stm-results-excel
description: Write an Excel rundown of every image in an STM analysis run, including interpretation, evidence, uncertainty, and failures, and save it in that run's overall-results folder. Use after analysis or when a run summary is requested.
---

# STM Results Excel

Create a readable `.xlsx` rundown for every input image in a specified STM run. Save it as `<run>/overall-results/image_results.xlsx`. This reporting skill does not itself authorize rerunning analysis or changing classifications.

## Read the run evidence

Read `run_manifest.json`, `overall-results/combined_analysis.json` or CSV, and each image's whole-image evidence JSON, sliced summary, and final evidence summary where available. Use the manifest's input list as the completeness baseline. Include failed, excluded, interrupted, and missing-result inputs, rather than silently dropping them. Match records by stable source identifier, not row order; keep same-stem images distinct.

Follow the available `stm-double-tip-analysis` interpretation principles and `saving-protocol` archive requirements. For workbook authoring, apply the available Spreadsheets skill and supported runtime. An existing `stm-double-tip-analysis/scripts/build_run_summary_workbook.mjs` may be adapted if available, but verify its columns and logic against this skill; it is not an independent source of conclusions.

## Workbook contents

Use two worksheets when group totals help: `Summary` and `Image Results`. Put the plain-language image rundown first among the detailed columns. Avoid one worksheet per image.

`Image Results` must have one row per input image with these core columns:

- Image filename or unique identifier.
- Analysis status: successful, failed, excluded, incomplete, or not analyzed.
- Reported overall interpretation and a concise, image-specific rundown.
- Whole-image classification and sliced/window classification, kept separate.
- Evidence supporting the result; limitations, conflicting evidence, or review needed.
- Relative paths or verified relative hyperlinks to the image evidence folder and final evidence summary.
- Error or exclusion reason, if applicable.

Write the rundown in 2–4 short sentences explaining what was found, what supports it, whether it is localized or scan-wide when known, and why uncertainty remains. Use concrete evidence from that image; avoid copying generic method descriptions into every row. `No clear evidence` does not establish a single tip. An unsuccessful analysis is not a negative finding.

Include useful quantitative fields when recorded: confidence (with its meaning), selected feature-map displacement dx/dy in pixels, physical displacement only with valid calibration, accepted pair count, independent spatial regions, ACF strength/prominence, pair confidence, window count, min/max/mean sliced score, transition window, partial-window indices, and scan/acquisition direction. Keep units in headers and numerical cells numeric. Do not label a raw score as a probability. Leave unavailable numerical values blank and explain missing evidence; never substitute zero.

Use the source's supported final classification without silently strengthening or merging it. Preserve `Likely double tip`, `Likely multiple tip`, `Possible/uncertain`, and `No clear evidence` when provided. If a source only reports `Likely double/multiple tip`, retain that wording and mark subtype unspecified. A sliced transition alone or a high ACF score does not establish a confirmed tip artifact. Flag disagreements between combined results and underlying evidence for review; distinguish any reviewer interpretation from the original automated output. Preserve user corrections separately with provenance.

`Summary` should identify the run and analysis/report timestamps separately, then show total inputs, successful/failed/excluded/incomplete counts, classification counts for analyzed images, and a brief group interpretation with limitations. Counts must reconcile to the input list. Do not rank the strongest image solely by maximum sliced score, or count missing/failed results as `No clear evidence`. If classification counts exclude unresolved records, show that exclusion explicitly.

## Formatting and traceability

Use a filterable detail table, frozen headers, wrapped narrative cells, readable widths/heights, and restrained header styling. Distinguish uncertain, failed, and incomplete rows without relying on color alone. Keep long evidence text readable rather than shrinking fonts to fit a single page. Add concise source paths and the relevant JSON field or record identifier where useful for tracing values. Use relative links from the workbook's `overall-results` location, normally `../images/...`, and verify their targets. Do not embed every image in the workbook; keep full-resolution evidence in each image folder.

Use formulas for summary counts when supported, recalculate, and check that formulas reflect classifications and analysis status correctly. Treat source text beginning with formula characters as literal text. Retain original classifications when correcting a reporting or formatting problem.

## Save and verify

Create `overall-results` if needed. Preserve an earlier workbook before replacing it; use `image_results - revision 2.xlsx` and subsequent revision numbers unless the user requests overwrite. For unfinished runs, clearly mark the workbook provisional and record its coverage timestamp; do not mark the analysis complete merely because the report exported.

Before delivery:

1. Reconcile image rows to the manifest and source list; resolve duplicates and omissions. Handle an empty run without invalid ranges or fabricated rows.
2. Spot-check a positive, uncertain, negative, and failed/missing image when those cases exist against underlying evidence. Confirm units, vector selection, missing values, and classification wording.
3. Verify summary totals, calculate formulas, and scan for unexpected formula errors.
4. Render and visually inspect both worksheets with the available spreadsheet tools. Fix clipped text, unreadable widths, and broken links; disclose unavailable checks rather than claiming they passed.
5. Export the actual `.xlsx` into this run's `overall-results`. Confirm it exists and opens with expected sheets and records where tooling permits.
6. Add its relative path to the run manifest's overall outputs without changing the analysis completion status. Preserve the report builder and normalized inputs under the run's materials or overall-results as appropriate, and refresh the saving-protocol inventory after the files and manifest stop changing.

Return the workbook location and briefly describe coverage and any unresolved source-evidence issues. Creating this skill is separate from executing it: generate a workbook only when the user's reporting request or analysis workflow calls for it.
