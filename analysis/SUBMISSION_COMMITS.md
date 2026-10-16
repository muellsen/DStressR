# Submission Commit Guide

Run these commands from the repository root. They deliberately exclude local
data, generated analysis results, internal manuscript audits, and unrelated
exploratory directories. No commits or push have been made by Codex.

```sh
cd /Users/cmueller/Documents/DStressR
git status --short
git diff --cached --stat
```

Start with an empty index; inspect any pre-existing staged changes before
continuing. Each group below should be reviewed before committing.

## 1. Inference and Defaults

Includes Welch--Satterthwaite centered uncertainty, single-pass variance
moderation, full-interaction centered-contrast covariance, total-effect
defaults, EVC parameter-conflict checks, single-pass moderation of fitted EVC
contrasts, regenerated R help, and regression tests.

```sh
git add -- R/model.R R/plots.R R/mixture.R man \
  tests/testthat/test-centered-inference.R tests/testthat/test-destress.R \
  tests/testthat/test-evc-adjustment.R
git diff --cached --check
git diff --cached --stat
git commit -m "Correct contrast uncertainty and default to total-effect inference"
```

## 2. Examples and Paper Reproducibility

Includes the two new total-effect preparers required by the canonical builder.
Specific-effect supplementary figures remain explicit and reproducible;
historical full-interaction specific-effect artifacts require regeneration
before reuse.

```sh
git add -- README.md \
  vignettes/DStressR.Rmd \
  vignettes/ecoli-promoter-screen.Rmd \
  vignettes/pfibox-promoter-activity-screen.Rmd \
  analysis/SUBMISSION_COMMITS.md \
  analysis/PAPER_FIGURES.md \
  analysis/build_paper_figures.R \
  analysis/ecoli_promoter_screen/INTERACTION_INFERENCE_NOTES.md \
  analysis/ecoli_promoter_screen/prepare_total_effect_discovery.R \
  analysis/ecoli_promoter_screen/prepare_interaction_total_effects.R \
  analysis/ecoli_promoter_screen/run_additive_evc_analysis.R \
  analysis/ecoli_promoter_screen/ADDITIVE_EVC_ANALYSIS.md \
  analysis/ecoli_promoter_screen/plot_discovery_figures.R \
  analysis/ecoli_promoter_screen/plot_interaction_figures.R \
  analysis/ecoli_promoter_screen/run_interaction_model.R \
  analysis/french_pfibox/plot_clustering_and_inference.R \
  analysis/french_pfibox/plot_french_cluster_alluvial.R \
  analysis/french_pfibox/plot_pfibox_model_cluster_and_discovery.R \
  analysis/french_pfibox/plot_pfibox_moderation_diagnostic.R \
  analysis/dryad_global_regulators/plot_rank_adjustment_figures.R \
  analysis/dryad_global_regulators/plot_weak_stress_three_target_inference.R
git diff --cached --check
git diff --cached --stat
git commit -m "Align public examples and paper figure workflows with total effects"
```

## 3. Manuscript

Includes the revised exposition, total-effect applications, supplementary
ordering, references, concise reporter--perturbation terminology, layout, and
explicit PFIbox clustering details. Rendered PDF/HTML remain local under the
existing manuscript policy.

```sh
git add -- paper/dstressr_package_manuscript/manuscript.qmd paper/references.bib
git diff --cached --check
git diff --cached --stat
git commit -m "Finalize total-effect manuscript and supplementary analyses"
```

## 4. Generated Documentation

Rebuild the local site to include the latest EVC documentation. GitHub Pages is
deployed by the existing pkgdown GitHub Action after pushing to a configured
branch; local rebuilding alone does not publish the site.

```sh
Rscript -e 'pkgdown::build_site(new_process = FALSE, install = TRUE, preview = FALSE)'
git add -- docs
git diff --cached --check
git diff --cached --stat
git commit -m "Rebuild package reference and tutorials for submission"
```

## Review and Push

```sh
git log -4 --oneline
git status --short
git push
```

Unrelated exploratory files may still appear as untracked. Do not add them
with `git add .`. Keep raw/private data and `analysis/outputs/` excluded.
