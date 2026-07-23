source(file.path("analysis", "_helpers.R"))

if (!requireNamespace("base64enc", quietly = TRUE)) {
  stop("Package `base64enc` is required to build the standalone HTML report.", call. = FALSE)
}

out_dir <- analysis_output_dir("binsfeld")
plot_script <- analysis_path("analysis", "ecoli_promoter_screen", "plot_comparison_figures.R")

required_files <- file.path(out_dir, c(
  "hit_overlap_venn.png",
  "pvalue_histograms.png",
  "effect_histograms.png",
  "effect_scatter.png",
  "pvalue_scatter_zoom.png",
  "binsfeld_destress_all_pair_comparison.tsv",
  "binsfeld_destress_significant_union.tsv",
  "destress_default_growth_exponents.tsv"
))
if (any(!file.exists(required_files))) {
  message("Missing one or more plot/table outputs; regenerating Binsfeld figures first.")
  status <- system2(file.path(R.home("bin"), "Rscript"), plot_script)
  if (!identical(status, 0L)) {
    stop("Could not regenerate Binsfeld plot outputs.", call. = FALSE)
  }
}

all_pairs <- read.delim(
  file.path(out_dir, "binsfeld_destress_all_pair_comparison.tsv"),
  check.names = FALSE,
  stringsAsFactors = FALSE
)
union_pairs <- read.delim(
  file.path(out_dir, "binsfeld_destress_significant_union.tsv"),
  check.names = FALSE,
  stringsAsFactors = FALSE
)
growth_exponents <- read.delim(
  file.path(out_dir, "destress_default_growth_exponents.tsv"),
  check.names = FALSE,
  stringsAsFactors = FALSE
)

fmt_num <- function(x, digits = 3) {
  ifelse(is.na(x), "", format(signif(as.numeric(x), digits), scientific = TRUE, trim = TRUE))
}
html_escape <- function(x) {
  x <- as.character(x)
  x <- gsub("&", "&amp;", x, fixed = TRUE)
  x <- gsub("<", "&lt;", x, fixed = TRUE)
  x <- gsub(">", "&gt;", x, fixed = TRUE)
  x <- gsub('"', "&quot;", x, fixed = TRUE)
  x
}
img_uri <- function(path) {
  paste0("data:image/png;base64,", base64enc::base64encode(path))
}
html_table <- function(df, columns = names(df), max_rows = Inf) {
  d <- df[, columns, drop = FALSE]
  if (is.finite(max_rows) && nrow(d) > max_rows) {
    d <- d[seq_len(max_rows), , drop = FALSE]
  }
  header <- paste0("<th>", html_escape(names(d)), "</th>", collapse = "")
  rows <- apply(d, 1, function(row) {
    paste0("<tr>", paste0("<td>", html_escape(row), "</td>", collapse = ""), "</tr>")
  })
  paste0("<table><thead><tr>", header, "</tr></thead><tbody>", paste(rows, collapse = "\n"), "</tbody></table>")
}

count_summary <- data.frame(
  Metric = c(
    "Promoter-compound pairs tested",
    "Binsfeld-style significant pairs",
    "DStressR significant pairs",
    "Overlap",
    "Binsfeld-only",
    "DStressR-only",
    "Neither"
  ),
  Count = c(
    nrow(all_pairs),
    sum(all_pairs$binsfeld_hit),
    sum(all_pairs$destress_hit),
    sum(all_pairs$overlap_class == "Both"),
    sum(all_pairs$overlap_class == "Binsfeld only"),
    sum(all_pairs$overlap_class == "DStressR only"),
    sum(all_pairs$overlap_class == "Neither")
  ),
  stringsAsFactors = FALSE
)
count_value <- function(metric) {
  count_summary$Count[match(metric, count_summary$Metric)]
}
binsfeld_hit_n <- count_value("Binsfeld-style significant pairs")
destress_hit_n <- count_value("DStressR significant pairs")
overlap_n <- count_value("Overlap")

promoter_summary <- aggregate(
  cbind(binsfeld_hit, destress_hit) ~ promoter,
  all_pairs,
  sum
)
names(promoter_summary) <- c("Promoter", "Binsfeld hits", "DStressR hits")
overlap_by_promoter <- aggregate(
  overlap_class ~ promoter,
  all_pairs,
  function(x) sum(x == "Both")
)
promoter_summary$Overlap <- overlap_by_promoter$overlap_class[
  match(promoter_summary$Promoter, overlap_by_promoter$promoter)
]
promoter_summary <- promoter_summary[order(promoter_summary$Promoter), ]

both <- union_pairs[union_pairs$overlap_class == "Both", ]
both$binsfeld_padj <- fmt_num(both$binsfeld_padj)
both$specific_padj_by_promoter <- fmt_num(both$specific_padj_by_promoter)
both$mean_z <- fmt_num(both$mean_z)
both$specific_effect <- fmt_num(both$specific_effect)

appendix <- union_pairs
appendix$mean_z <- fmt_num(appendix$mean_z)
appendix$binsfeld_pvalue <- fmt_num(appendix$binsfeld_pvalue)
appendix$binsfeld_padj <- fmt_num(appendix$binsfeld_padj)
appendix$specific_effect <- fmt_num(appendix$specific_effect)
appendix$specific_pvalue <- fmt_num(appendix$specific_pvalue)
appendix$specific_padj_by_promoter <- fmt_num(appendix$specific_padj_by_promoter)

binsfeld_only <- appendix[appendix$overlap_class == "Binsfeld only", ]
destress_only <- appendix[appendix$overlap_class == "DStressR only", ]
overlap_hits <- appendix[appendix$overlap_class == "Both", ]

growth_table <- growth_exponents
for (col in intersect(c("a_raw", "a_raw_se", "alpha_raw", "alpha_raw_se", "alpha_shrunk"), names(growth_table))) {
  growth_table[[col]] <- fmt_num(growth_table[[col]])
}

report_file <- file.path(out_dir, "binsfeld_destress_shareable_report.html")
created <- format(Sys.time(), "%Y-%m-%d %H:%M %Z")

html <- paste0(
'<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<title>DStressR comparison on the Binsfeld et al. reporter screen</title>
<style>
body { font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif; line-height: 1.48; color: #111827; margin: 0; background: #f8fafc; }
main { max-width: 1060px; margin: 0 auto; padding: 34px 28px 56px; background: white; }
h1 { font-size: 32px; margin: 0 0 8px; }
h2 { margin-top: 34px; border-top: 1px solid #e5e7eb; padding-top: 24px; }
h3 { margin-top: 24px; }
p, li { font-size: 15px; }
.meta { color: #4b5563; margin-bottom: 26px; }
.callout { background: #eff6ff; border-left: 4px solid #2563eb; padding: 12px 16px; margin: 18px 0; }
.grid { display: grid; grid-template-columns: repeat(2, minmax(0, 1fr)); gap: 18px; }
.figure { margin: 22px 0; }
.figure img { width: 100%; border: 1px solid #e5e7eb; background: white; }
.caption { color: #4b5563; font-size: 13px; margin-top: 6px; }
table { border-collapse: collapse; width: 100%; margin: 12px 0 22px; font-size: 13px; }
th, td { border: 1px solid #e5e7eb; padding: 6px 8px; vertical-align: top; }
th { background: #f3f4f6; text-align: left; }
code { background: #f3f4f6; padding: 1px 4px; border-radius: 3px; }
.small { font-size: 13px; color: #4b5563; }
@media print { body { background: white; } main { max-width: none; padding: 20px; } .grid { grid-template-columns: 1fr; } }
</style>
</head>
<body><main>
<h1>DStressR comparison on the Binsfeld et al. reporter screen</h1>
<div class="meta">Self-contained local report generated ', html_escape(created), ' from the DStressR repository.</div>

<div class="callout">
<strong>Headline result.</strong> The reproduced Binsfeld-style WT analysis calls ', binsfeld_hit_n, ' promoter-compound hits. The default DStressR modeled-response analysis calls ', destress_hit_n, ' hits. The exact promoter-compound overlap is ', overlap_n, ' hits.
</div>

<h2>Data and Analysis Scope</h2>
<p>This report uses the public <em>E. coli</em> reporter-screen data from Binsfeld et al. The DStressR package ships an AUC-level table, <code>binsfeld_reporter_auc</code>, prepared from the PLOS S3 Data supplement, and an author score/Z-score table, <code>binsfeld_reporter_scores</code>, prepared from the PLOS S4 Data supplement.</p>
<p>The comparison is limited to WT reporter rows. For the DStressR analysis, rows marked <code>removed == "No"</code> are used, water controls are collapsed to <code>Water</code>, promoter-specific growth-response exponents are estimated from control wells, and the empty-vector reporter is excluded from the default modeled-response testing set.</p>

<h2>Methods Compared</h2>
<h3>Reproduced Binsfeld-style rule</h3>
<ul>
<li>Input: WT Z-scores from <code>binsfeld_reporter_scores</code>.</li>
<li>For each promoter and compound, compare replicate/concentration Z-scores against water controls using a Wilcoxon test.</li>
<li>Adjust p-values within promoter using Benjamini-Hochberg.</li>
<li>Call a hit when promoter-wise adjusted p-value &lt; 0.05 and absolute mean Z-score &gt; 1.</li>
</ul>
<h3>DStressR model</h3>
<ul>
<li>Input: WT AUC rows from <code>binsfeld_reporter_auc</code> with <code>removed == "No"</code>.</li>
<li>Response: package-default modeled response, <code>log2(lux_auc) - alpha_g * log2(od_auc)</code>, with promoter-specific <code>alpha_g</code> estimated from water controls.</li>
<li>Model: DStressR model preset, empirical-Bayes moderation enabled, technical terms for replicate and dose level.</li>
<li>Call a hit using promoter-specific effect p-values adjusted within promoter at FDR 0.05.</li>
</ul>

<h3>Estimated Growth-Response Exponents</h3>
<p>The table below records the default response-modeling step used for this screen.</p>',
html_table(
  growth_table,
  columns = intersect(c("promoter", "alpha_raw", "alpha_raw_se", "alpha_shrunk", "n_controls", "shrink_weight", "alpha_covariates"), names(growth_table))
),
'
<h2>Numerical Summary</h2>',
html_table(count_summary),
'<h3>Hit Counts by Promoter</h3>',
html_table(promoter_summary),

'<h2>Figures</h2>
<div class="figure"><img alt="Hit overlap Venn diagram" src="', img_uri(file.path(out_dir, "hit_overlap_venn.png")), '"><div class="caption">Figure 1. Exact promoter-compound hit overlap.</div></div>
<div class="figure"><img alt="P-value histograms" src="', img_uri(file.path(out_dir, "pvalue_histograms.png")), '"><div class="caption">Figure 2. Raw and promoter-wise BH-adjusted p-value distributions. Colored bars are method-specific hits.</div></div>
<div class="figure"><img alt="Effect histograms" src="', img_uri(file.path(out_dir, "effect_histograms.png")), '"><div class="caption">Figure 3. Effect-score distributions for original mean Z-scores and DStressR specific effects.</div></div>
<div class="grid">
<div class="figure"><img alt="Effect scatter" src="', img_uri(file.path(out_dir, "effect_scatter.png")), '"><div class="caption">Figure 4. Effect-score comparison across all tested promoter-compound pairs.</div></div>
<div class="figure"><img alt="P-value scatter zoomed" src="', img_uri(file.path(out_dir, "pvalue_scatter_zoom.png")), '"><div class="caption">Figure 5. Zoomed raw p-value comparison.</div></div>
</div>

<h2>Interpretation</h2>
<p>The two analyses agree on a shared core of ', overlap_n, ' promoter-compound interactions. DStressR calls reflect modeled growth-adjusted AUC responses, promoter-specific effects after compound-wide centering, and moderated model-based uncertainty rather than the original Z-score/Wilcoxon rule. Binsfeld-only hits tend to be pairs whose mean Z-score and Wilcoxon evidence pass the original thresholds but whose DStressR promoter-specific effect is weaker after model adjustment.</p>

<h2>Top Overlapping Hits</h2>',
html_table(
  both[order(as.numeric(both$binsfeld_padj)), ],
  columns = c("promoter", "compound", "mean_z", "binsfeld_padj", "specific_effect", "specific_padj_by_promoter", "destress_hit_class"),
  max_rows = 20
),

'<h2>Appendix: Union of Significant Pairings</h2>
<p>The method-specific unique hit tables are written as TSV files alongside this report:</p>
<ul>
<li><code>analysis/outputs/binsfeld/binsfeld_only_significant_pairs.tsv</code></li>
<li><code>analysis/outputs/binsfeld/destress_only_significant_pairs.tsv</code></li>
<li><code>analysis/outputs/binsfeld/overlapping_significant_pairs.tsv</code></li>
</ul>

<h3>Binsfeld-only Significant Pairs</h3>',
html_table(
  binsfeld_only,
  columns = c("promoter", "compound", "mean_z", "binsfeld_padj", "binsfeld_direction", "specific_effect", "specific_padj_by_promoter")
),

'<h3>DStressR-only Significant Pairs</h3>',
html_table(
  destress_only,
  columns = c("promoter", "compound", "mean_z", "binsfeld_padj", "specific_effect", "specific_padj_by_promoter", "destress_hit_class")
),

'<h3>Overlapping Significant Pairs</h3>',
html_table(
  overlap_hits,
  columns = c("promoter", "compound", "mean_z", "binsfeld_padj", "binsfeld_direction", "specific_effect", "specific_padj_by_promoter", "destress_hit_class")
),

'<h3>Full Union of Significant Pairings</h3>
<p>The table below is the union of all promoter-compound pairs significant by either method. The same table is written as a shareable TSV at <code>analysis/outputs/binsfeld/binsfeld_destress_significant_union.tsv</code>.</p>',
html_table(
  appendix,
  columns = c("promoter", "compound", "overlap_class", "mean_z", "binsfeld_padj", "binsfeld_direction", "specific_effect", "specific_padj_by_promoter", "destress_hit_class")
),

'<h2>Reproducibility</h2>
<p>Regenerate figures and TSVs from the repository root with:</p>
<pre><code>Rscript analysis/ecoli_promoter_screen/plot_comparison_figures.R
Rscript analysis/ecoli_promoter_screen/build_shareable_report.R</code></pre>
<p class="small">Primary public sources: PLOS Biology article DOI 10.1371/journal.pbio.3003260 and Zenodo DOI 10.5281/zenodo.15600688.</p>
</main></body></html>'
)

writeLines(html, report_file, useBytes = TRUE)
message("Wrote standalone report: ", report_file)
