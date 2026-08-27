# Zaslaver et al. promoter activity and OD time courses

Long-form version of Dataset S1 from Zaslaver et al. (2009). The source
file contains, for each growth condition, a promoter-activity matrix
followed by the matching OD matrix. These values are author-processed
promoter activity and OD summaries, not the raw fluorescence
plate-reader traces.

## Usage

``` r
zaslaver_promoter_timecourse
```

## Format

A data frame with 622,080 rows and 9 columns:

- condition:

  Machine-readable condition label.

- condition_label:

  Source condition label.

- reporter_index:

  Source row index of the reporter construct.

- reporter_id:

  Unique reporter identifier derived from row index and source row
  index.

- reporter:

  Original source reporter label.

- time_index:

  Source measurement-column index.

- time_min:

  Nominal time in minutes, assuming 16-minute sampling intervals as
  reported by Zaslaver et al.

- promoter_activity:

  Author-processed promoter activity.

- od:

  Optical-density value from the matching OD block.

## Source

Zaslaver et al. (2009), PLOS Computational Biology, Dataset S1.
