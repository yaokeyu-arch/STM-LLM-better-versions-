import fs from "node:fs/promises";
import { SpreadsheetFile, Workbook } from "@oai/artifact-tool";

const inputPath = process.argv[2];
const outputPath = process.argv[3];
if (!inputPath || !outputPath) throw new Error("Usage: node build_run_summary_workbook.mjs input.json output.xlsx");

const rows = JSON.parse(await fs.readFile(inputPath, "utf8"));
const workbook = Workbook.create();
const summary = workbook.worksheets.add("Summary");
const detail = workbook.worksheets.add("Image Results");
const font = "Arial";

summary.showGridLines = false;
summary.getRange("A2:F2").merge();
summary.getRange("A2").values = [["STM double/multiple-tip analysis"]];
summary.getRange("A2").format.font = { name: font, size: 15, bold: true, color: "#17365D" };
summary.getRange("A3:F3").format.borders = { bottom: { style: "thin", color: "#4472C4" } };
summary.getRange("A5:B9").values = [
  ["Measure", "Result"],
  ["Images analyzed", rows.length],
  ["Likely double/multiple tip", null],
  ["Possible/uncertain", null],
  ["No clear evidence", null],
];
summary.getRange("B7").formulas = [[`=COUNTIF('Image Results'!$D$2:$D$${rows.length + 1},"Likely double/multiple tip")`]];
summary.getRange("B8").formulas = [[`=COUNTIF('Image Results'!$D$2:$D$${rows.length + 1},"Possible/uncertain")`]];
summary.getRange("B9").formulas = [[`=COUNTIF('Image Results'!$D$2:$D$${rows.length + 1},"No clear evidence")`]];
summary.getRange("A5:B5").format = { fill: "#17365D", font: { name: font, bold: true, color: "#FFFFFF" }, horizontalAlignment: "center" };
summary.getRange("A6:B9").format.font = { name: font, size: 10 };
summary.getRange("A11:F11").values = [["Interpretation", "Primary evidence", "Strongest image", "Maximum sliced score", "Run folder", "Status"]];
const strongest = [...rows].sort((a, b) => (b.max_score ?? -1) - (a.max_score ?? -1))[0];
summary.getRange("A12:F12").values = [[
  "Classification combines whole-image morphology with sliced evidence; autocorrelation alone is insufficient.",
  "Symmetric ACF candidates, duplicated morphology, spatial/sliced consistency, height weighting, and lattice rejection.",
  strongest?.file ?? "",
  strongest?.max_score ?? null,
  outputPath.split("/overall-results/")[0],
  rows.some((r) => r.error) ? "Completed with errors" : "Complete",
]];
summary.getRange("A11:F11").format = { fill: "#4472C4", font: { name: font, bold: true, color: "#FFFFFF" }, horizontalAlignment: "center", wrapText: true };
summary.getRange("A12:F12").format = { font: { name: font, size: 10 }, wrapText: true, verticalAlignment: "top" };
summary.getRange("D12").format.numberFormat = "0.000";
summary.getRange("A1:F12").format.font.name = font;
summary.getRange("A1:F12").format.autofitColumns();
summary.getRange("A1:F12").format.autofitRows();
summary.getRange("A12:B12").format.columnWidth = 48;
summary.getRange("C:C").format.columnWidth = 22;
summary.getRange("D:D").format.columnWidth = 22;
summary.getRange("E12").format.columnWidth = 46;
summary.getRange("F:F").format.columnWidth = 16;
summary.getRange("A12:F12").format.rowHeight = 58;

const headers = ["Image", "Whole-image classification", "Sliced classification", "Overall interpretation", "Reason", "Scan direction", "Acquisition order", "Windows", "Minimum score", "Maximum score", "Mean score", "Final-window score", "Change window", "Partial windows", "Image evidence folder", "Error"];
const values = rows.map((r) => [
  r.file, r.whole_classification, r.sliced_classification, r.combined_classification,
  r.whole_reason, r.scan_direction, r.acquisition_order, r.window_count,
  r.min_score, r.max_score, r.mean_score, r.last_score,
  r.change_window_index, Array.isArray(r.partial_correlation_window_indices) ? r.partial_correlation_window_indices.join(", ") : "",
  r.processed_npz ? r.processed_npz.replace(/\/[^/]+$/, "") : "", r.error,
]);
detail.getRangeByIndexes(0, 0, 1, headers.length).values = [headers];
detail.getRangeByIndexes(1, 0, values.length, headers.length).values = values;
detail.showGridLines = false;
detail.freezePanes.freezeRows(1);
detail.getRangeByIndexes(0, 0, 1, headers.length).format = { fill: "#17365D", font: { name: font, bold: true, color: "#FFFFFF" }, horizontalAlignment: "center", verticalAlignment: "center", wrapText: true };
detail.getRangeByIndexes(1, 0, values.length, headers.length).format = { font: { name: font, size: 10 }, verticalAlignment: "top" };
detail.getRange(`I2:L${rows.length + 1}`).format.numberFormat = "0.000";
detail.getRange(`A1:P${rows.length + 1}`).format.autofitColumns();
detail.getRange(`A1:P${rows.length + 1}`).format.autofitRows();
detail.getRange(`B2:E${rows.length + 1}`).format.wrapText = true;
detail.getRange("E:E").format.columnWidth = 42;
detail.getRange("O:O").format.columnWidth = 48;
detail.getRange(`D2:D${rows.length + 1}`).conditionalFormats.add("containsText", { text: "Likely", format: { fill: "#E2F0D9", font: { color: "#006100", bold: true } } });
detail.getRange(`D2:D${rows.length + 1}`).conditionalFormats.add("containsText", { text: "No clear", format: { fill: "#F2F2F2", font: { color: "#666666" } } });

workbook.recalculate();
console.log((await workbook.inspect({ kind: "table", sheetId: "Summary", range: "A2:F12", include: "values,formulas", tableMaxRows: 15, tableMaxCols: 8 })).ndjson);
console.log((await workbook.inspect({ kind: "match", searchTerm: "#REF!|#DIV/0!|#VALUE!|#NAME\\?|#N/A|#NUM!|#NULL!|#SPILL!|#CALC!", options: { useRegex: true, maxResults: 100 }, summary: "formula error scan" })).ndjson);
const preview = await workbook.render({ sheetName: "Summary", range: "A1:F12", scale: 1.5, format: "png" });
await fs.writeFile(outputPath.replace(/\.xlsx$/i, "_preview.png"), new Uint8Array(await preview.arrayBuffer()));
await fs.mkdir(outputPath.replace(/\/[^/]+$/, ""), { recursive: true });
const output = await SpreadsheetFile.exportXlsx(workbook);
await output.save(outputPath);
