# Binsfeld et al. reporter screen AUC data

A public *E. coli* reporter-screen data set from Binsfeld et al. (2025),
prepared as an AUC-level long table for DStressR examples and tests. The
rows are reporter/strain/replicate/well observations from the PLOS
Biology S3 Data supplement. `drug` keeps the original Binsfeld
perturbation labels, including `Water_1` and `Water_2`. `compound` is a
derived grouping column that collapses the water control wells to
`Water` for summaries. `dose_level` is derived from
`concentration_index` so that larger values correspond to higher
compound concentration.

## Usage

``` r
binsfeld_reporter_auc
```

## Format

A data frame with 24,576 rows and 13 columns:

- strain:

  Reporter host strain.

- reporter:

  Reporter label, including `EVC`.

- replicate:

  Reporter replicate number.

- well:

  384-well plate coordinate.

- drug:

  Original drug/control label from the source table.

- compound:

  Derived drug/control grouping label, with water controls collapsed.

- concentration_index:

  Dose-series index from the source table.

- dose_level:

  Dose-oriented serial dilution level; larger values correspond to
  higher concentration.

- concentration_ug_ml:

  Compound concentration in micrograms per ml.

- od_auc:

  Optical-density area under the curve.

- lux_auc:

  Luminescence area under the curve.

- od_auc_per_lux_auc:

  Source-table OD/LUX AUC ratio.

- removed:

  Author quality-control flag.

## Source

Binsfeld et al. (2025), PLOS Biology, S3 Data.

## Details

The source article is https://doi.org/10.1371/journal.pbio.3003260. The
associated Zenodo code/data archive is
https://doi.org/10.5281/zenodo.15600688.
