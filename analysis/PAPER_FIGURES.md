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
| 5; total-effect promoter p-value supplement | `ecoli_promoter_screen/plot_discovery_figures.R --target=total` | `binsfeld_total_effect_discovery/binsfeld_discovery.rds` |
| Reporter-relative effect discovery supplement | `ecoli_promoter_screen/plot_discovery_figures.R` | `binsfeld_evc_calibrated/binsfeld_discovery.rds` |
| 6 | `french_pfibox/build_selected_activity_figure.R --target=total` | Integrated-activity results, clustering, and saved 20-pair trajectories |
| PFIbox reporter-relative-effect supplement | `french_pfibox/plot_pfibox_model_cluster_and_discovery.R` | Same standard result table and saved clustering |
| Reporter-relative effect matrices | `ecoli_promoter_screen/plot_binsfeld_effect_matrix_views.R` | `binsfeld_effect_paper.rds` |
| Interaction-model total-effect supplement | `ecoli_promoter_screen/plot_interaction_figures.R` | `ecoli_interaction_total_effects/binsfeld_interaction_paper.rds` |
| PFIbox total-effect moderation supplement | `french_pfibox/plot_pfibox_moderation_diagnostic.R --target=total` | Standard ordinary/moderated results TSVs |
| Kochanowski supplement | `kochanowski_promoter_scaling/plot_paper_figures.R` | `kochanowski_paper.rds` |
| Dash effect decomposition | `dryad_global_regulators/plot_rank_adjustment_figures.R` | `dash_rank_paper.rds` |
| Dash inference supplement | `dryad_global_regulators/plot_weak_stress_three_target_inference.R` | Standard default/rank-1 results TSVs |

Figure 6 assembles the total-effect heatmap, four manually selected activity
trajectories, moderated p-value histogram, and four-labelled-pair volcano.
Its underlying renderers are `plot_pfibox_model_cluster_and_discovery.R` and
`plot_clustering_and_inference.R`. Its diagnostic producer is
`french_pfibox/compute_paper_diagnostics.R`; its annotation producer is
`map_french_timecourse_clusters_to_destress.R`. All use the same explicit run
tag and cluster suffix. Files are never chosen by modification time. Published
French cluster assignments are annotations, not re-estimated time-course clusters.

The figure-only contract is checked by
`Rscript analysis/check_paper_figure_contract.R`. It rejects direct fitting,
inference, decomposition, or clustering calls in the retained plotting
scripts. Summaries used only for display (histogram bins, labels, hit overlaps,
color limits) remain in the figure generators. Paper-specific panel assembly
is preserved; reusable statistical computations remain in DStressR.
Label placement uses a fixed seed. Optional non-manuscript matrix views,
including on-the-fly clustered heatmaps, are available separately through
`ecoli_promoter_screen/plot_exploratory_matrix_views.R` and are not part of the
canonical figure-only pipeline.

## Total-Effect Discovery Comparison

The total-effect version of Figure 5 reuses the saved Binsfeld fits:

```sh
Rscript analysis/ecoli_promoter_screen/prepare_total_effect_discovery.R
Rscript analysis/ecoli_promoter_screen/plot_discovery_figures.R --target=total
```

Outputs are isolated in `analysis/outputs/binsfeld_total_effect_discovery/`.
The original-workflow results are unchanged; both DStressR analyses use
moderated total-effect tests and within-reporter BH adjustment at 0.05.
The preparer records the input fit paths and checks the stored BH values.
The default `--target=relative` retains the supplementary reporter-relative-effect figure
and its output directory. The paper builder prepares the total-effect artifact
in its diagnostics stage and renders both targets in its figure stage.
Neither command edits the manuscript or refits the models. The preparer also
exports separate total-effect union and candidate-annotation tables; the
existing reporter-relative-effect lists remain unchanged.

## Interaction-Model Total Effects

Supplementary Note F compares the same total-effect target across reporter-wise
and full interaction models:

```sh
Rscript analysis/ecoli_promoter_screen/prepare_interaction_total_effects.R
Rscript analysis/ecoli_promoter_screen/plot_interaction_figures.R
```

The preparer reads the saved full interaction-model result artifact, including
its raw total-effect variances, and applies the package's empirical-Bayes
routine with the full model's residual degrees of freedom. Earlier artifacts
stored ordinary total-effect tests; their p-values are checked against the
saved standard errors before reconstructing moderated tests. No response or
coefficient is refitted. Newly generated interaction fits are also saved as RDS
files by `run_interaction_model.R`.

Outputs go to `analysis/outputs/ecoli_interaction_total_effects/`. The earlier
reporter-relative-effect outputs remain in `analysis/outputs/ecoli_interaction_model/`.
Regenerate them with the corrected joint-covariance implementation before
rendering them with `--target=relative`. The paper builder
prepares the total comparison in its diagnostics stage. A local snapshot before
the change is stored in
`analysis/outputs/inference_audit_2026-10-16/pre_interaction_total_effects.tar.gz`.

## PFIbox Total-Effect Analysis

Figure 6 uses total effects of the signed integral of reconstructed unscaled
activity. Run the annotated stages in `french_pfibox/README.md`, then:

```sh
Rscript analysis/french_pfibox/build_selected_activity_figure.R --target=total
```

Outputs are isolated in
`analysis/outputs/french_pfibox/temporal_unscaled_integral_analysis/selected_four_figure6/`.
Panel a uses correlation distance and complete linkage on both axes. Panel b
shows ompC / Ampicillin-18 and msrB / Sulfamethoxazole-18 on the left, with
recN / Ciprofloxacin-18 and rmf / NalidixicAcid-18 on the right. These are
manually chosen from the documented joint-score 10-positive / 10-negative
display selection. Panels c and d use `total_pvalue`, `total_padj_global` and
`total_effect`; the same four pairs are marked in the heatmap and volcano.
The manuscript example has 59,301 tests and 5,101 global BH calls at 0.05.

Supplementary Note H uses ordinary/moderated total-effect diagnostics and
reporter-relative inference from this same integral fit. The paper builder
sets the explicit source directory, run tag, global adjustment scope and
clustering suffix. The figure stage does not fit models or recompute clusters.
Abandoned PFIbox response-summary, PCA and alternative-clustering explorations
have been removed; the exact PA reconstruction audit remains available.

## Preserved Baseline

The pre-refactor sources and principal figure outputs are stored locally in
`analysis/outputs/inference_audit_2026-10-16/pre_figure_refactor.tar.gz`.
The earlier pre-inference-revision archives remain unchanged. Archives and
generated artifacts are ignored by Git; source generators and this contract
are versioned.
