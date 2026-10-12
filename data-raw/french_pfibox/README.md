# French PFIbox Data Provenance

This directory is reserved for local source and processed files from the PFIbox
promoter-reporter screen of French, Coutts, and Brown (2018), *Cell Systems*,
doi:10.1016/j.cels.2018.07.004.

## Source

The primary deposited dataset is listed as Mendeley Data
doi:10.17632/fmv6vv4bsf.1. The public code repository is
<https://github.com/sfrench007/pfibox>, which contains the PFIbox scripts and
sample ImageJ output files used to validate parsers.

The local DStressR workflow expects the full multi-treatment `gfpData.Rdata`
object from the deposited data, plus optional supplementary cluster annotation
tables from the article. These files are external source data and are not
committed to Git.

## Local Layout

Place downloaded source files under:

```text
data-raw/french_pfibox/source/
```

The current local analysis used:

```text
data-raw/french_pfibox/source/gfpData.Rdata
data-raw/french_pfibox/source/french_pfibox_table_s3_clusterdata_mmc2.xlsx
```

Derived audit tables and compact intermediate files are written to:

```text
data-raw/french_pfibox/processed/
```

The `source/` and `processed/` directories are ignored by Git because they can
contain large raw or derived files. Re-run the scripts in
`analysis/french_pfibox/` to regenerate processed tables from the downloaded
source files.
