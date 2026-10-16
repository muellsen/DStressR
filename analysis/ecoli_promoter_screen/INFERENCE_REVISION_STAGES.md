# Inference revision: 2026-10-16

The pre-change submission is retained under
`analysis/outputs/inference_audit_2026-10-16/`: two archives preserve sources,
manuscript, fitted objects, and all four application output directories.
Do not overwrite these archives.

- [x] Audit candidate code on unchanged fits; compare all four applications.
- [x] Promote and document single moderation, effective df, and unequal-df prior moments.
- [x] Align main statistical framework and Supplementary Note C.
- [x] Regenerate Binsfeld results and paper figures; verify unchanged effects.
- [x] Regenerate PFIbox results, ordinary-test comparison, and paper figures.
- [x] Regenerate Kochanowski results/figures and correct summary dimensions.
- [x] Regenerate Dash results/figures, including the archived rank-test df discrepancy.
- [x] Run package tests, render manuscript, and check numerical/text consistency.

Verification: `Rscript analysis/ecoli_promoter_screen/verify_inference_promotion.R`
checks 41 regenerated quantities against the independent audit. Source and
figure checksums are saved beside the protected archives. PDF and HTML renders
complete successfully. The regression tests include unequal-df prior moments,
Welch-Satterthwaite centering, and once-only moderation.

## Result Changes

| Application | Reporter-wise FDR calls after revision |
| --- | --- |
| Binsfeld, no EVC / EVC | 80 / 92 (unchanged) |
| PFIbox, specific | 5,184 (previously 5,103); global 4,123 (previously 3,988) |
| Kochanowski | Unchanged, exploratory summary-level scores |
| Dash, total / specific / rank-adjusted total | 36 / 42 / 30 |

The Dash archive had 26 rank-adjusted calls because its p-values used
unaugmented residual df. Both pre-change source code and the promoted code
give 30; this separate historical discrepancy is documented in the audit.
All fitted effect estimates remain unchanged.

No automatic Git commit or push is part of this workflow.
