# Reporter-Reporter-relative Inference Audit

The audited statistical update has been promoted to the package, manuscript,
and application outputs. Progress is recorded in `INFERENCE_REVISION_STAGES.md`.
The pre-change submission remains available in the archives below.

## Preserved Baseline

`analysis/outputs/inference_audit_2026-10-16/prechange_submission_snapshot.tar.gz`
preserves the package sources, tests, manuscript QMD/PDF, references, Binsfeld
drivers, and principal Binsfeld output directories from before the change.
`baseline_R/` contains the frozen source functions, and `baseline_fits.rds`
contains both fitted Binsfeld objects. These are local, ignored artifacts.

Do not recreate this checkpoint from changed sources. A clean checkout needs
the checkpoint restored before running the historical comparison.

## Promoted Calculation

The reporter-wise path constructs centered-effect uncertainty from raw total
variances, computes Welch--Satterthwaite effective degrees of freedom, and
moderates the centered variance once. The empirical prior uses per-contrast
log-chi-square corrections for unequal degrees of freedom. Total-effect
moderation remains a separate calculation. Result tables expose raw variances
and raw/final reporter-relative-effect degrees of freedom for audit.

The full interaction-model path retains its common residual degrees of freedom.
The dedicated screen-specific workflows and supplied-summary path are unchanged.
The separate exported `global_effect` diagnostic still uses a representative
median residual df. It is not the specific- or total-effect test used for the
manuscript's hit calls.

## Commands

From the repository root:

```sh
Rscript analysis/ecoli_promoter_screen/audit_single_moderation.R
Rscript analysis/ecoli_promoter_screen/check_single_moderation_null.R
Rscript analysis/ecoli_promoter_screen/audit_other_examples_moderation.R
Rscript analysis/ecoli_promoter_screen/verify_inference_promotion.R
quarto render paper/dstressr_package_manuscript/inference_revision.qmd --to pdf
Rscript -e 'devtools::test(reporter="summary")'
```

All comparison outputs go to `analysis/outputs/inference_audit_2026-10-16/`.
The null check defaults to 1,000 simulated screens per scenario; set
`DSTRESSR_NULL_SIMULATIONS` to change this count.
Historical comparisons read their original tables from the protected archives,
so regenerating current application outputs does not change the baseline.

## Findings

- All 658 fitted effects in each Binsfeld variant remain identical.
- The candidate retains the same 80/92 within-reporter FDR calls without/with EVC.
- The archived no-EVC p-values used the earlier minimum-df convention; an
  explicit reconstruction recovers them to numerical precision. The prechange
  source uses per-reporter df and has a small numerical difference without any
  change in hit calls. The archived EVC results match the prechange source.
- Gaussian null checks show broadly nominal raw rejection rates in large-df
  designs. The small heterogeneous design retains modest excess rejection under
  candidate moderation, so it should not be described as an exact calibrated test.
- The focused PDF documents equations, numerical comparisons, and limitations.
  Application outputs are regenerated with the promoted calculation; the
  archived outputs retain the historical results.

## Other Examples

The three additional output directories are preserved in
`other_examples_prechange_outputs.tar.gz`. Derived inputs and results are
checksum-checked before and after the comparison. Cached fits and all candidate
outputs are under the audit's `other_examples/` directory.

- PFIbox specific calls: 5,103 to 5,184 reporter-wise (172 gained, 91 lost),
  and 3,988 to 4,123 globally (168 gained, 33 lost). Total-effect calls and
  fitted matrices/clusters are unchanged.
- Dash total/specific calls remain 36/42. Rank-adjusted total calls are 30
  in both the source baseline and candidate, versus 26 in the archived figure.
  The archived p-values are exactly reproduced using unaugmented residual df;
  this is a pre-existing difference independent of the candidate change.
- Kochanowski is numerically unchanged. Its saved input actually has 26
  conditions including glucose, so its 64 reporters give 1,600 non-reference
  contrasts; these dimensions have been corrected in the manuscript.
