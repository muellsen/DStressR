# French PFIbox Promoter-Reporter Analysis

This folder contains a DStressR-oriented audit and analysis scaffold for the
PFIbox dataset of French, Coutts, and Brown (2018), *Cell Systems*,
doi:10.1016/j.cels.2018.07.004.

The assay differs from the Binsfeld luminescence screen. The public workflow
contains colony-level GFP fluorescence time courses and promoter-position
metadata for the Zaslaver *E. coli* promoter-GFP library. Growth is represented
through a constitutive-reporter biomass proxy, not through matched OD600
kinetics for every reporter colony. The primary DStressR workflow in this folder
therefore treats the deposited promoter activity as a supplied response and
summarizes it over a user-defined time window.

## Data Sources

Primary article data:

- Mendeley Data doi:10.17632/fmv6vv4bsf.1

Code and sample files:

- <https://github.com/sfrench007/pfibox>

The GitHub sample contains one untreated time-course file and is useful for
parser validation. A multi-treatment `gfpData.Rdata` object or equivalent
processed time-course output is needed for DStressR perturbation analysis.

## Scripts

Validate the public GitHub sample parser:

```sh
Rscript analysis/french_pfibox/prepare_pfibox_sample_data.R
```

Audit the full deposited `gfpData.Rdata` object:

```sh
Rscript analysis/french_pfibox/audit_pfibox_data.R
```

Run the DStressR alpha-comparison analysis after placing a multi-treatment
`gfpData.Rdata` file under `data-raw/french_pfibox/source/`:

```sh
Rscript analysis/french_pfibox/run_french_pfibox_destress_analysis.R
```

The script pools the deposited `*-0` treatments as the reference condition,
keeps nonzero drug-dose entries as perturbations, and runs the model twice:
once with reporter-specific response exponents estimated from the reference
colonies and once with the fixed exponent `alpha = 1`.
By default, it reuses the cached compact AUC table when available. Set
`PFIBOX_REBUILD_AUC=true` to rebuild that table from `gfpData.Rdata`.
For a quick smoke test, set `PFIBOX_REPORTER_LIMIT` and
`PFIBOX_PERTURBATION_LIMIT`; subset outputs are tagged so they do not overwrite
the full-run files.

Run the promoter-activity workflow on a reproducible subset:

```sh
PFIBOX_REPORTER_LIMIT=250 \
PFIBOX_PERTURBATION_LIMIT=8 \
PFIBOX_T_START_H=2 \
PFIBOX_T_END_H=18 \
Rscript analysis/french_pfibox/run_pfibox_promoter_activity_subset.R
```

This script uses
`y_i^{PA}(T) = |T|^{-1} int_T PA_i(t) dt` as the supplied DStressR response,
where the time window `T` is controlled by `PFIBOX_T_START_H` and
`PFIBOX_T_END_H`. The defaults analyze roughly 2,000 reporter-perturbation
pairs before missing-value filtering and write the response table, pair-level
DStressR results, perturbation diagnostics, and a compact summary figure.

Generate paper-style plots for the same promoter-activity subset:

```sh
PFIBOX_RUN_TAG=reporters250_perturbations8_T2-18h \
Rscript analysis/french_pfibox/plot_pfibox_promoter_activity_figures.R
```

These figures mirror the Binsfeld application figures after response
construction: total and specific effect matrices, a rank-variance diagnostic,
specific-effect p-value histogram, volcano plot, and a table-style display of
the strongest specific effects. The response-modeling parameter panels are not
generated for this workflow because promoter activity is supplied directly.

Generated compact tables are written under
`data-raw/french_pfibox/processed/`; figures and analysis outputs are written
under `analysis/outputs/french_pfibox/`. Both locations are ignored by Git.
