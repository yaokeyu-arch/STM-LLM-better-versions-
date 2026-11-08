import fs from "node:fs/promises";
import { FileBlob, SpreadsheetFile, Workbook } from "@oai/artifact-tool";

const root = "/Users/cindy/Documents/Codex/2026-09-22/se";
const outputDir = process.env.OUTPUT_DIR || `${root}/outputs/stm_evaluation`;
const referencePath = "/Users/cindy/Desktop/topo_desc.xlsx";
const blindPath = process.env.ANALYSIS_PATH || `${root}/outputs/blind_analysis_double_tip_images/combined_analysis.json`;
const filePrefix = process.env.FILE_PREFIX || "stm";
const font = "Arial";

const blind = JSON.parse(await fs.readFile(blindPath, "utf8"));
const referenceBook = await SpreadsheetFile.importXlsx(await FileBlob.load(referencePath));
const referenceRows = referenceBook.worksheets.getItem("Sheet1").getRange("A1:E855").values;
const reference = new Map(referenceRows.slice(1).map((r) => [String(r[0]), {
  file: String(r[0]), tear: Number(r[1] || 0), noise: Number(r[2] || 0),
  lattice: Number(r[3] || 0), doubleTip: Number(r[4] || 0),
}]));

function parseCandidate(text = "") {
  const re = /^\(([+-]?\d+),([+-]?\d+)\) px\s+\(([+-]?[\d.]+),([+-]?[\d.]+)\) nm\s+ACF F\/B=([\d.-]+)\/([\d.-]+)\s+prominence=([\d.-]+) sigma\s+corr F\/B=([\d.-]+)\/([\d.-]+)\s+feature-F1 F\/B=([\d.-]+)\/([\d.-]+)\s+correlated tiles F\/B=(\d+)\/16\/(\d+)\/16\s+feature tiles F\/B=(\d+)\/16\/(\d+)\/16(?:\s+periodicity penalty=([\d.-]+))?$/;
  const m = text.match(re);
  if (!m) return {};
  const n = m.slice(1).map(Number);
  return {
    dxPx:n[0], dyPx:n[1], dxNm:n[2], dyNm:n[3], acfF:n[4], acfB:n[5],
    prominence:n[6], corrF:n[7], corrB:n[8], f1F:n[9], f1B:n[10],
    corrTilesF:n[11], corrTilesB:n[12], featureTilesF:n[13], featureTilesB:n[14],
  };
}

function blindClass(row) {
  if (String(row.combined_classification).startsWith("Likely")) return "Likely double/multiple tip";
  if (String(row.whole_classification).startsWith("Likely")) return "Likely double/multiple tip";
  if (String(row.sliced_classification).includes("Partial double-tip correlations") ||
      String(row.sliced_classification).includes("Persistent double-tip evidence from scan start") ||
      String(row.sliced_classification).includes("Persistent single-to-double-tip transition")) {
    return "Likely double/multiple tip";
  }
  if (String(row.whole_classification).startsWith("Possible") ||
      String(row.sliced_classification).includes("single-to-double-tip transition")) {
    return "Possible/uncertain";
  }
  return "No clear evidence";
}

const compared = blind.map((b) => {
  const ref = reference.get(b.file);
  if (!ref) throw new Error(`No reference row for ${b.file}`);
  const prediction = blindClass(b);
  const predictedPositive = prediction === "Likely double/multiple tip";
  let outcome;
  if (ref.doubleTip === 1 && predictedPositive) outcome = "Correct positive";
  else if (ref.doubleTip === 0 && prediction === "No clear evidence") outcome = "Correct negative";
  else if (ref.doubleTip === 0 && predictedPositive) outcome = "False positive";
  else if (ref.doubleTip === 1 && prediction === "Possible/uncertain") outcome = "Inconclusive against positive reference";
  else if (ref.doubleTip === 1) outcome = "False negative";
  else outcome = "Inconclusive against negative reference";
  return {...b, ...ref, prediction, outcome, candidate: parseCandidate(b.whole_top_candidate)};
});

const headers = [
  "File", "Blind classification", "Reference double tip", "Outcome",
  "Tear", "Noise", "Lattice", "Δx (px)", "Δy (px)", "Δx (nm)", "Δy (nm)",
  "ACF forward", "ACF backward", "Peak prominence (σ)",
  "Morphology correlation forward", "Morphology correlation backward",
  "Irregular-feature F1 forward", "Irregular-feature F1 backward",
  "Correlated regions forward", "Correlated regions backward",
  "Feature-matched regions forward", "Feature-matched regions backward",
  "Sliced result", "Maximum sliced score", "Explanation",
];

function record(row) {
  const c = row.candidate;
  const explanation = row.outcome === "Correct positive"
    ? "Blind analysis and reference both identify double-tip imaging."
    : row.outcome.startsWith("Inconclusive")
      ? "Blind evidence was suggestive but did not meet the threshold for a confident positive classification."
      : row.outcome === "False negative"
        ? "Reference identifies double tip, but morphology-confirmed evidence did not meet the blind-analysis threshold."
        : row.outcome === "False positive"
          ? "Blind analysis identified double tip, but the reference label is negative."
          : "Blind analysis and reference both indicate no double-tip classification.";
  return [
    row.file, row.prediction, row.doubleTip ? "Yes" : "No", row.outcome,
    row.tear, row.noise, row.lattice,
    c.dxPx ?? null, c.dyPx ?? null, c.dxNm ?? null, c.dyNm ?? null,
    c.acfF ?? null, c.acfB ?? null, c.prominence ?? null,
    c.corrF ?? null, c.corrB ?? null, c.f1F ?? null, c.f1B ?? null,
    c.corrTilesF ?? null, c.corrTilesB ?? null, c.featureTilesF ?? null, c.featureTilesB ?? null,
    row.sliced_classification, Number(row.max_score), explanation,
  ];
}

function colName(index) {
  let n = index + 1, s = "";
  while (n) { n--; s = String.fromCharCode(65 + (n % 26)) + s; n = Math.floor(n / 26); }
  return s;
}

async function createWorkbook(title, subtitle, rows, filename, tabColor) {
  const wb = Workbook.create();
  const sheet = wb.worksheets.add("Results");
  sheet.showGridLines = false;
  sheet.tabColor = tabColor;
  const lastCol = colName(headers.length - 1);
  sheet.getRange(`A2:${lastCol}2`).merge();
  sheet.getRange("A2").values = [[title]];
  sheet.getRange(`A3:${lastCol}3`).merge();
  sheet.getRange("A3").values = [[subtitle]];
  sheet.getRange(`A5:${lastCol}5`).values = [headers];
  if (rows.length) sheet.getRange(`A6:${lastCol}${5 + rows.length}`).values = rows.map(record);
  else {
    sheet.getRange(`A6:${lastCol}6`).merge();
    sheet.getRange("A6").values = [["No records in this category."]];
  }

  const usedEnd = Math.max(6, 5 + rows.length);
  sheet.getRange(`A2:${lastCol}${usedEnd}`).format.font = {name:font, size:10, color:"#1F2937"};
  sheet.getRange(`A2:${lastCol}2`).format.font = {name:font, size:15, bold:true, color:"#16324F"};
  sheet.getRange(`A3:${lastCol}3`).format.font = {name:font, size:10, italic:true, color:"#4B5563"};
  sheet.getRange(`A5:${lastCol}5`).format = {
    fill:"#1F4E78", font:{name:font,size:10,bold:true,color:"#FFFFFF"},
    horizontalAlignment:"center", verticalAlignment:"center", wrapText:true,
    borders:{preset:"inside",style:"thin",color:"#FFFFFF"},
  };
  sheet.getRange(`A6:${lastCol}${usedEnd}`).format.verticalAlignment = "center";
  if (rows.length) {
    sheet.getRange(`A6:${lastCol}${usedEnd}`).format.borders = {
      insideHorizontal:{style:"thin",color:"#D7DEE7"},
      bottom:{style:"thin",color:"#D7DEE7"},
    };
    sheet.getRange(`E6:G${usedEnd}`).format.numberFormat = "0";
    sheet.getRange(`H6:I${usedEnd}`).format.numberFormat = "0";
    sheet.getRange(`J6:K${usedEnd}`).format.numberFormat = "0.000";
    sheet.getRange(`L6:R${usedEnd}`).format.numberFormat = "0.000";
    sheet.getRange(`S6:V${usedEnd}`).format.numberFormat = "0";
    sheet.getRange(`X6:X${usedEnd}`).format.numberFormat = "0.000";
    sheet.getRange(`A6:A${usedEnd}`).format.font = {name:font,size:10,bold:true,color:"#16324F"};
    sheet.getRange(`Y6:Y${usedEnd}`).format.wrapText = true;
    sheet.getRange(`D6:D${usedEnd}`).conditionalFormats.add("containsText", {
      text:"False", format:{fill:"#FDE2E2",font:{color:"#9B1C1C",bold:true}},
    });
    sheet.getRange(`D6:D${usedEnd}`).conditionalFormats.add("containsText", {
      text:"Correct", format:{fill:"#E4F3E8",font:{color:"#166534",bold:true}},
    });
    sheet.getRange(`D6:D${usedEnd}`).conditionalFormats.add("containsText", {
      text:"Inconclusive", format:{fill:"#FFF4CE",font:{color:"#8A5A00",bold:true}},
    });
    sheet.tables.add(`A5:${lastCol}${usedEnd}`, true, `${filename.replace(/[^A-Za-z0-9]/g, "")}Table`);
  }
  sheet.freezePanes.freezeRows(5);
  sheet.freezePanes.freezeColumns(1);
  sheet.getRange(`A5:${lastCol}${usedEnd}`).format.autofitColumns();
  sheet.getRange(`A5:${lastCol}${usedEnd}`).format.autofitRows();
  const widths = {
    A:150,B:185,C:105,D:215,E:55,F:55,G:55,H:65,I:65,J:70,K:70,L:88,M:88,N:105,
    O:115,P:115,Q:115,R:115,S:100,T:100,U:115,V:115,W:220,X:105,Y:370,
  };
  for (const [col,width] of Object.entries(widths)) sheet.getRange(`${col}:${col}`).format.columnWidthPx = width;
  sheet.getRange("2:2").format.rowHeightPx = 26;
  sheet.getRange("3:3").format.rowHeightPx = 22;
  sheet.getRange("5:5").format.rowHeightPx = 58;
  if (rows.length) sheet.getRange(`6:${usedEnd}`).format.rowHeightPx = 54;
  wb.recalculate();

  const inspect = await wb.inspect({kind:"table",sheetId:"Results",range:`A2:Y${usedEnd}`,include:"values,formulas",tableMaxRows:25,tableMaxCols:25,maxChars:16000});
  console.log(inspect.ndjson);
  const errors = await wb.inspect({kind:"match",searchTerm:"#REF!|#DIV/0!|#VALUE!|#NAME\\?|#N/A|#NUM!|#NULL!|#SPILL!|#CALC!",options:{useRegex:true,maxResults:100},summary:`${filename} formula error scan`});
  console.log(errors.ndjson);
  const preview = await wb.render({sheetName:"Results",range:`A1:Y${Math.min(usedEnd,12)}`,scale:1,format:"png"});
  await fs.writeFile(`${outputDir}/${filename}.png`,new Uint8Array(await preview.arrayBuffer()));
  const out = await SpreadsheetFile.exportXlsx(wb);
  await out.save(`${outputDir}/${filename}.xlsx`);
}

await fs.mkdir(outputDir,{recursive:true});
const correct = compared.filter((r) => r.outcome.startsWith("Correct"));
const falsePositives = compared.filter((r) => r.outcome === "False positive");
const falseNegatives = compared.filter((r) => r.outcome === "False negative" || r.outcome.startsWith("Inconclusive against positive"));

await createWorkbook("Correct STM classifications", `${correct.length} predictions agree with the reference label.`, correct, `${filePrefix}_correct_classifications`, "#2F6B4F");
await createWorkbook("STM false positives", `${falsePositives.length} scans were called double tip when the reference label was negative.`, falsePositives, `${filePrefix}_false_positives`, "#9B1C1C");
await createWorkbook("STM false negatives and inconclusive positives", `${falseNegatives.length} reference-positive scans were not confidently classified as double tip.`, falseNegatives, `${filePrefix}_false_negatives`, "#B26A00");

console.log(JSON.stringify({correct:correct.length,falsePositives:falsePositives.length,falseNegatives:falseNegatives.length,outputDir},null,2));
