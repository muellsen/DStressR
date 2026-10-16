# Paper figure pipeline

Run from the repository root, with the downloaded data described in each
dataset README available locally:

```sh
Rscript analysis/build_paper_figures.R --analysis
Rscript analysis/build_paper_figures.R --diagnostics
Rscript analysis/build_paper_figures.R --figures
quarto render paper/dstressr_package_manuscript/manuscript.qmd
```

`--all` runs all three stages. The default is `--figures`: it reads saved
artifacts, does not fit models or recalculate inference, and fails if an input
is missing. Individual generators can also be run directly.

## Artifact Contract

`_paper_artifacts.R` defines a versioned local RDS container with:

- `results`: named, unrenamed DStressR `results()` tables;
- `assays`: prepared response inputs where needed for response plots;
- `parameters`: response/calibration or rank parameters;
- `diagnostics`: saved matrix decompositions, clustering, and trend diagnostics;
- `annotations`: published scores/cluster labels, display orders, and metadata;
- `provenance`: run identifiers and diagnostic settings.

This container belongs to paper reproducibility, not to the installed package
API. It retains numerical data, not serialized plot objects or environments.
The existing fitted-object RDS files remain available separately.

## Figure Map

| Figure | Generator under `analysis/` | Primary saved input |
| --- | --- | --- |
| 1, 2 | Hand-edited illustration assets under `paper/figures/` | PNG/PDF assets |
| 3; response-distribution supplement | `ecoli_promoter_screen/plot_response_modeling.R` | `binsfeld_response_paper.rds` |
| 4 | `ecoli_promoter_screen/plot_binsfeld_effect_matrix_and_variance_diagnostic.R` | `binsfeld_effect_paper.rds` |
| 5; promoter p-value supplement | `ecoli_promoter_screen/plot_discovery_figures.R` | `binsfeld_discovery.rds` |
| 6 | `french_pfibox/plot_clustering_and_inference.R` | Explicit PFIbox diagnostic and mapping RDS files |
| Specific-effect matrices | `ecoli_promoter_screen/plot_binsfeld_effect_matrix_views.R` | `binsfeld_effect_paper.rds` |
| Interaction-model supplement | `ecoli_promoter_screen/plot_interaction_figures.R` | `binsfeld_interaction_paper.rds` |
| PFIbox moderation supplement | `french_pfibox/plot_pfibox_moderation_diagnostic.R` | Standard ordinary/moderated results TSVs |
| Kochanowski supplement | `kochanowski_promoter_scaling/plot_paper_figures.R` | `kochanowski_paper.rds` |
| Dash effect decomposition | `dryad_global_regulators/plot_rank_adjustment_figures.R` | `dash_rank_paper.rds` |
| Dash inference supplement | `dryad_global_regulators/plot_weak_stress_three_target_inference.R` | Standard default/rank-1 results TSVs |

Figure 6 assembles `plot_pfibox_model_cluster_and_discovery.R` and
`plot_french_cluster_alluvial.R`. Its diagnostic producer is
`french_pfibox/compute_paper_diagnostics.R`; its annotation producer is
`map_french_timecourse_clusters_to_destress.R`. All use the same explicit run
tag and cluster suffix. Files are never chosen by modification time. Published
French cluster assignments are annotations, not re-estimated time-course clusters.

The figure-only contract is checked by
`Rscript analysis/check_paper_figure_contract.R`. It rejects direct fitting,
inference, decomposition, or clustering calls in the 13 canonical plotting
scripts. Summaries used only for display (histogram bins, labels, hit overlaps,
color limits) remain in the figure generators. Paper-specific panel assembly
is preserved; reusable statistical computations remain in DStressR.
Label placement uses a fixed seed. Optional non-manuscript matrix views,
including on-the-fly clustered heatmaps, are available separately through
`ecoli_promoter_screen/plot_exploratory_matrix_views.R` and are not part of the
canonical figure-only pipeline.

## Preserved Baseline

The pre-refactor sources and principal figure outputs are stored locally in
`analysis/outputs/inference_audit_2026-10-16/pre_figure_refactor.tar.gz`.
The earlier pre-inference-revision archives remain unchanged. Archives and
generated artifacts are ignored by Git; source generators and this contract
are versioned.
