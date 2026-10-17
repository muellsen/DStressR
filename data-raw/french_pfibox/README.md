# French PFIbox Data Provenance

This directory is reserved for local source and processed files from the PFIbox
promoter-reporter screen of French, Coutts, and Brown (2018), *Cell Systems*,
doi:10.1016/j.cels.2018.07.004.

## Source

The primary deposited dataset is Mendeley Data version 1 (24 August 2018),
DOI <https://doi.org/10.17632/fmv6vv4bsf.1>, available at
<https://data.mendeley.com/datasets/fmv6vv4bsf/1>. Its licence is CC BY 4.0;
derived data exports must retain attribution and describe their transformations.
The public code repository is
<https://github.com/sfrench007/pfibox>, which contains the PFIbox scripts and
sample ImageJ output files used to validate parsers.

The local DStressR workflow expects the full multi-treatment `gfpData.Rdata`
object from the deposited data, plus the published Table S3 cluster annotations
for the supplementary comparison. These files are external source data and are not
committed to Git.

## Local Layout

Download `gfpData.Rdata` from the version-1 dataset page and place it at the
path below. Alternatively, run these commands from the repository root:

```sh
mkdir -p data-raw/french_pfibox/source
curl --fail --location --retry 3 \
  --output data-raw/french_pfibox/source/gfpData.Rdata \
  https://data.mendeley.com/public-files/datasets/fmv6vv4bsf/files/d7526382-74f1-485c-8678-ccb4d21f33f2/file_downloaded
shasum -a 256 data-raw/french_pfibox/source/gfpData.Rdata
```

Expected size: **1,311,202,733 bytes**. Expected SHA-256:

```text
41b3f09bc7a8c3c435af9d03d0e33cede018ba7a5fd0ec900cf6f7d07079ed63
```

The download URL and checksum were resolved from Mendeley's version-1 file
metadata. Use the dataset page if the direct endpoint changes; do not silently
substitute a different version or the GitHub protocol sample.

For the published temporal-cluster comparison, obtain Table S3 (`mmc2.xlsx`)
from the article's supplementary files and save it under the second path
below. The worksheet used is `ClusterData`.

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

The source object has 46 condition entries plus a common colony layout:
13 untreated arrays and 33 non-reference antibiotic-dose treatments. The
retained analysis reconstructs unscaled activity from `constit.normalized`;
the deposited `promoter.activity` matrices are already min-max-scaled.
The source-stage audit and exact reconstruction checks are documented in
`analysis/french_pfibox/FRENCH_FIGURE3_AUDIT.md`.
