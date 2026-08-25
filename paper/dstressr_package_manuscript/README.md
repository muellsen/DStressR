# DStressR Package Manuscript

This manuscript is the methodological and software paper for DStressR.

It supports the package narrative while remaining separate from the CRAN
package source. The package build ignores `paper/`; the manuscript is rendered
from this folder and reads figures/tables produced by scripts under
`analysis/`.

The current application section uses two public examples:

- the Binsfeld et al. *E. coli* promoter-compound screen, reproduced by
  `analysis/ecoli_promoter_screen/`;
- the Dash et al. fluorescence/OD regulator-stress experiment, reproduced by
  `analysis/dryad_global_regulators/`.

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

## Rebuild Workflow

Run commands from the repository root. The most important paper-generation
scripts are:

```sh
Rscript analysis/ecoli_promoter_screen/plot_modeling_step_figures.R
Rscript analysis/ecoli_promoter_screen/run_evc_calibrated_analysis.R
Rscript analysis/ecoli_promoter_screen/plot_binsfeld_effect_matrix_and_variance_diagnostic.R
Rscript analysis/ecoli_promoter_screen/run_interaction_sensitivity.R
Rscript analysis/dryad_global_regulators/run_weak_stress_window_pipeline.R
Rscript analysis/dryad_global_regulators/run_weak_stress_window_calibrated_pipeline.R
Rscript analysis/dryad_global_regulators/run_weak_stress_window_rank1_adjustment.R
Rscript analysis/dryad_global_regulators/plot_weak_stress_three_target_inference.R
quarto render paper/dstressr_package_manuscript/manuscript.qmd
```

Generated outputs under `analysis/outputs/` are intentionally ignored by Git.
They can be regenerated from the scripts above when the package code changes.
