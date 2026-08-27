# Zaslaver et al. reporter annotation classes

Reporter annotation class matrix prepared from Dataset S3 of Zaslaver et
al. (2009). Annotation columns are binary indicators from the source
workbook. As in
[zaslaver_promoter_activity](https://muellsen.github.io/DStressR/reference/zaslaver_promoter_activity.md),
`reporter_id` is the unique reporter construct identifier and `reporter`
preserves the source label.

## Usage

``` r
zaslaver_promoter_annotations
```

## Format

A data frame with 1,920 rows and 249 columns. The first columns are
`reporter_index`, `reporter_id`, and `reporter`; remaining columns are
binary annotation indicators.

## Source

Zaslaver et al. (2009), PLOS Computational Biology, Dataset S3.
