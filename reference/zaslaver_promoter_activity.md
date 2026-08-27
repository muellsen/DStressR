# Zaslaver et al. fixed-growth-rate promoter activity data

Public *E. coli* promoter-activity data from Zaslaver et al. (2009),
prepared from Dataset S2 of the PLOS Computational Biology article. Each
row is one reporter/condition/growth-rate summary. The `reporter` column
keeps the original source label, while `reporter_id` combines the source
row index and label because several reporter labels occur more than once
in the source workbooks and two reporter labels vary across condition
sheets.

## Usage

``` r
zaslaver_promoter_activity
```

## Format

A data frame with 23,040 rows and 10 columns:

- reporter_index:

  Source row index of the reporter construct.

- reporter_id:

  Unique reporter identifier derived from row index and source row
  index.

- reporter:

  Original source reporter label.

- growth_rate:

  Fixed growth rate at which promoter activity was summarized, in
  divisions per hour.

- condition:

  Machine-readable condition label.

- condition_label:

  Source condition label.

- promoter_activity:

  Promoter activity for the reporter under the condition at the fixed
  growth rate.

- mean_promoter_activity:

  Mean promoter activity across conditions at the fixed growth rate.

- sd_promoter_activity:

  Standard deviation of promoter activity across conditions at the fixed
  growth rate.

- cv_promoter_activity:

  Coefficient of variation of promoter activity across conditions at the
  fixed growth rate.

## Source

Zaslaver et al. (2009), PLOS Computational Biology, Dataset S2.

## Details

The source article is https://doi.org/10.1371/journal.pcbi.1000545. The
supplementary Excel files were obtained from the CaltechAUTHORS record
https://authors.library.caltech.edu/records/g76ja-xga32.
