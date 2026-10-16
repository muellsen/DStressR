# DStressR Package Manuscript

This is the statistical-methods manuscript for DStressR.

It supports the package narrative while remaining separate from the CRAN
package source. The package build ignores `paper/`; the manuscript is rendered
from this folder and reads figures/tables produced by scripts under
`analysis/`.

The current application section uses two public examples:

- the Binsfeld et al. *E. coli* promoter-compound screen, reproduced by
  `analysis/ecoli_promoter_screen/`;
- the French et al. PFIbox fluorescence screen, reproduced by
  `analysis/french_pfibox/`.

Supplementary Note G contains the Kochanowski promoter-activity compendium
(`analysis/kochanowski_promoter_scaling/`) and the Dash et al. fluorescence/OD
regulator-stress experiment (`analysis/dryad_global_regulators/`).

Campylobacter and other exploratory workflows are not part of this manuscript
track unless they are explicitly promoted in a later revision.

## Render

From the repository root:

```sh
quarto render paper/dstressr_package_manuscript/manuscript.qmd
```

## Primary Inputs

- Package code under `R/`.
- Package-shipped public datasets:
  - `binsfeld_reporter_auc`
  - `binsfeld_reporter_scores`
- Public *E. coli* analysis scripts under `analysis/ecoli_promoter_screen/`.
- Public fluorescence/OD regulator-stress scripts under
  `analysis/dryad_global_regulators/`.
- Shared bibliography in `paper/references.bib`.
- Hand-assembled manuscript figures in `paper/figures/`.
- Generated figure and table outputs under `analysis/outputs/`.

## Build Workflow

Run these stages from the repository root:

```sh
Rscript analysis/build_paper_figures.R --analysis
Rscript analysis/build_paper_figures.R --diagnostics
Rscript analysis/build_paper_figures.R --figures
quarto render paper/dstressr_package_manuscript/manuscript.qmd
```

The first stage prepares responses, fits models, and saves standardized results
and model parameters. The second computes diagnostics and attaches published
annotations. The third only renders figures from explicit saved artifacts.
Use `--all` for all stages. The complete figure map and artifact schema are
documented in `analysis/PAPER_FIGURES.md`.

Generated outputs under `analysis/outputs/` are intentionally ignored by Git.
They can be regenerated from the scripts above when the package code changes.

The data provenance and download instructions live in each analysis directory.
PFIbox model-only reruns can use `PFIBOX_RESPONSE_TABLE` with a response TSV
already produced for the same subset and window. The model drivers save fit
objects separately from inference tables so `results()` can be rerun directly.

The October 2026 inference revision and protected pre-change submission are
documented in `analysis/ecoli_promoter_screen/INFERENCE_REVISION_STAGES.md`
and `README_inference_audit.md`. The historical comparison requires the local
archives; a clean checkout can reproduce the current analysis without them.
