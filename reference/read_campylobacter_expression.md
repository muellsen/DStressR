# Read exported Campylobacter expression data

Read exported Campylobacter expression data

## Usage

``` r
read_campylobacter_expression(expression_file, libmap_file)
```

## Arguments

- expression_file:

  Path to `expression_values.tsv.gz`.

- libmap_file:

  Path to `LibMap.txt`.

## Value

A data frame with `srn_code` and `ProductName` joined in.
