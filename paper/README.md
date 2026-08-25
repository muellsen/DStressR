# Manuscripts

This directory separates manuscript assembly from the R package and from
dataset-specific analysis scripts. The CRAN package does not include this
folder; `.Rbuildignore` excludes `paper/`.

The current canonical manuscript is the DStressR package/methods paper in
`dstressr_package_manuscript/`. Other manuscript drafts are retained locally
unless explicitly promoted.

## Files

- `dstressr_package_manuscript/`: DStressR method/package manuscript. This is
  the primary paper and uses the public *E. coli* promoter-compound screen plus
  the fluorescence/OD regulator-stress example as applications.
- `references.bib`: shared bibliography for manuscript rendering.
- `figures/`: hand-assembled or manually edited manuscript figures. Paper
  figures generated from data should instead be produced by scripts under
  `analysis/` and written to `analysis/outputs/`.
- `campylobacter_analysis_manuscript/`: earlier Campylobacter analysis draft,
  retained as local/future manuscript material.
- `manuscript.qmd` and `MANUSCRIPT_SPLIT_PLAN.md`: legacy combined draft and
  split notes, retained locally for provenance.

## Render

From the repository root:

```sh
quarto render paper/dstressr_package_manuscript/manuscript.qmd
```

The main manuscript depends on generated figures and tables under
`analysis/outputs/`. Rebuild these from the dataset-specific workflows before
rendering when outputs are absent or stale.

The earlier Campylobacter and combined drafts may still render locally, but
they are not part of the canonical DStressR package-manuscript workflow.
