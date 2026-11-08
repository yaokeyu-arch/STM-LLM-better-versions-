---
name: saving-protocol
description: Preserve all input images, working materials, intermediate evidence, and final outputs for STM analysis runs, reruns, and exports in a self-contained, non-overwriting run archive.
---

# Saving Protocol

Apply this skill whenever running, rerunning, or exporting STM analysis. Archive the actual material used, not only the final figures or links to source files. This skill governs saving; it does not itself request an analysis run.

## Create a run

Use `scripts/create_run_layout.py <source-folder> <workspace>/runs --images <input-filenames...>` to create a group-named run without overwriting an earlier run. The helper creates the basic layout and manifest; create the additional archive folders below as needed. Repeated group names receive ` - run 2`, ` - run 3`, etc.; generic names use `run 1`, `run 2`, etc. Honor an explicitly selected destination.

Keep a separate folder for every input image, including failed or excluded images. Preserve relative source paths when stems collide, and record that mapping. Never modify source data. Copy originals byte-for-byte into the run before processing where possible; retain the source path, byte size, and SHA-256 checksum. If archival copying fails, record the problem and do not claim the archive is complete.

## Required archive

```text
<run name>/
  run_manifest.json
  file_inventory.json
  materials/
    code/                 # exact scripts and local modules actually used
    skills/               # applicable skill instructions at run time
    references/           # supplied reference images, masks, annotations, calibration, labels
    environment/          # runtime and relevant dependency versions
    commands.jsonl        # commands, working directories, settings, timestamps, exit codes
  overall-results/
    combined_analysis.json
    combined_analysis.csv
    # group comparisons, spreadsheets, plots, interpretation
  images/
    <image identifier>/
      source/             # original SXM/image files and associated metadata
      # original render, flattened and noise-reduced images, noise mask
      # ACF data and annotated maps; feature pairing and morphology comparisons
      # sliced/window data and maps; numerical evidence JSON/CSV/NPZ
      # final evidence summary; stdout/stderr and failure logs
```

Keep all outputs and calculations associated with one image inside its image folder, including intermediate and rejected-candidate evidence produced by the process. Put only group-level outputs in `overall-results/`. Preserve native-resolution figures and underlying numerical data; a montage or screenshot does not replace them. Retain native editable outputs where generated, such as workbooks and plotting code.

Copy all supporting material actually used: user-provided annotations and corrections, parameter/configuration files, calibration data, relevant reference documents, and exact analysis scripts with their local dependencies. Record their role and original location. Do not indiscriminately copy unrelated directories, credentials, environment secrets, caches, or virtual environments. Record relevant software versions instead of archiving installed packages. Preserve provenance for downloaded material when used.

Record effective parameters, random seeds where applicable, channels, scan direction, preprocessing choices, input-to-output relationships, and any manual decisions. Log execution commands and errors without secrets. Write logs and manifest progress during the run so interruption does not erase the record. When retrying, preserve earlier outputs and logs in separate attempt folders or start a new run; never silently replace them.

Ground-truth labels must stay unread until blind predictions are saved when the analysis skill requires blind evaluation. Archive labels afterward for evaluation, with their role and access stage recorded. Preserve user-verified corrections and distinguish them from blind predictions.

## Inventory and completion

Create `file_inventory.json` after outputs stop changing. Include every archived file except the inventory itself, with its run-relative path, byte size, SHA-256, role, and source or generating step where known. Hash files in chunks. Do not claim provenance that cannot be recovered. Use run-relative paths in manifests and output references so the archive remains usable if moved; record external source paths only as provenance.

Update `run_manifest.json` with:

- creation/completion timestamps and status (`running`, `complete`, `complete_with_failures`, or `interrupted`);
- input count, successes, failures, exclusions, and reasons;
- per-image source copies, evidence folders, and outcome;
- relative paths to overall results, inventory, code, configuration, logs, and environment records;
- any missing material or limitations affecting reproducibility.

Before reporting completion, verify every input is represented, original-copy checksums match, every inventoried file exists and its checksum matches, expected per-image evidence is present, and all result references resolve inside this run. Record failures in combined JSON/CSV and keep their logs. A successful script exit alone does not establish archive completeness. Mark incomplete archival checks explicitly and retain partial evidence.

For an existing run created before this protocol, backfill only material that is still available. Label later copies and reconstructed settings as such; do not assert they are the exact versions originally used without evidence.

Give the user the run location and links to overall results and the inventory. If a downloadable bundle is useful, package the verified archive without replacing the folder; preserve its directory structure. Follow the host chat's deliverable-location conventions.
