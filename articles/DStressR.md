# Get started with DStressR

`DStressR` models promoter-activity responses in high-throughput
chemical genomics screens. The package is designed as the
stress-response counterpart to
[`DGrowthR`](https://bio-datascience.github.io/DGrowthR/): DGrowthR
handles growth-curve modeling, while DStressR handles promoter-compound
effects after accounting for growth and technical structure.

This vignette uses a small simulated screen so that the model-based
named workflow can be run without private data files. DStressR also
exposes compatibility workflows through the same entry point:
`workflow = "median_polish"` for the original median-polish p-value
workflow and `workflow = "empty_vector_control"` for the
empty-vector-control workflow.

## Input files and table shape

DStressR expects a long measurement table with one row per measured
promoter-compound-replicate observation. In a complete rectangular
screen, the row count is approximately:

``` text
n_promoters x (n_compound_wells + n_control_wells) x n_replicates
```

The table may contain additional rows when the same screen is repeated
across batches, library plates, measurement plates, or experimental
days. The important point is that these technical covariates remain in
the table so they can be modeled rather than accidentally averaged away.

The required biological and technical information is:

- promoter or reporter construct
- compound or library-well identifier
- negative-control compound, usually `DMSO`
- luminescence summary
- growth summary
- replicate and optional technical covariates such as batch or plate

Column names are not fixed. They are mapped explicitly in
[`prepare_assay()`](https://muellsen.github.io/DStressR/reference/prepare_assay.md).

For the original Campylobacter promoter-library workflow, the
convenience helper
[`read_campylobacter_expression()`](https://muellsen.github.io/DStressR/reference/read_campylobacter_expression.md)
joins two exported files:

- `expression_values.tsv.gz`: one row per
  promoter-library-well-replicate observation. It must contain either
  `srn_code` or both `libplate` and `well`, plus the promoter,
  luminescence, growth, and technical columns used downstream.
- `LibMap.txt`: one row per compound/library well with columns
  `Library plate`, `Well`, `ProductName`, and `Catalog Number`.

The helper reconstructs the library-well key as
`paste0("lp", Library plate, "_", Well)` for `LibMap.txt`, joins the
compound annotations onto the expression table, and returns a long table
with the same number of rows as `expression_values.tsv.gz`.

``` r

expression_df <- read_campylobacter_expression(
  expression_file = "expression_values.tsv.gz",
  libmap_file = "LibMap.txt"
)

assay <- prepare_assay(
  expression_df,
  promoter = "promoter",
  compound = "srn_code",
  control = "DMSO",
  lux = "LUX.AUC_16",
  growth = "od_16h.measured",
  batch = "batch",
  plate = "libplate",
  replicate = "replicate"
)
```

## Simulate a chemical-genomics screen

``` r

library(DStressR)

screen <- simulate_screen(
  n_promoters = 12,
  n_compounds = 40,
  n_replicates = 3,
  seed = 1
)

head(screen)
#>   promoter compound replicate batch od_16h.measured LUX.AUC_16 truth_specific
#> 1       P1     DMSO        r1    b2       0.4116687   293.9048              0
#> 2       P2     DMSO        r1    b2       0.3906898   432.5738              0
#> 3       P3     DMSO        r1    b2       0.3333026   235.1525              0
#> 4       P4     DMSO        r1    b2       0.3000495   850.4902              0
#> 5       P5     DMSO        r1    b2       0.3410560   346.1222              0
#> 6       P6     DMSO        r1    b2       0.5140298   378.6924              0
#>   truth_global
#> 1            0
#> 2            0
#> 3            0
#> 4            0
#> 5            0
#> 6            0
```

The simulated table contains one row per promoter-compound-replicate
well. It includes luminescence, growth, batch, replicate, and hidden
truth columns used only to check the example.

``` r

names(screen)
#> [1] "promoter"        "compound"        "replicate"       "batch"          
#> [5] "od_16h.measured" "LUX.AUC_16"      "truth_specific"  "truth_global"
```

## Prepare the assay

[`prepare_assay()`](https://muellsen.github.io/DStressR/reference/prepare_assay.md)
converts luminescence and growth summaries into a growth-adjusted log2
response. By default, DStressR estimates promoter-specific growth
exponents from DMSO control wells and shrinks them toward a global
control-well slope.

``` r

assay <- prepare_assay(
  screen,
  promoter = "promoter",
  compound = "compound",
  control = "DMSO",
  lux = "LUX.AUC_16",
  growth = "od_16h.measured",
  batch = "batch",
  replicate = "replicate"
)

attr(assay, "destress")$growth_exponent_fit
#>    promoter control_n log_growth_sd a_raw a_raw_se a_raw_df alpha_raw
#> 1        P1         3     0.2789856    NA       NA       NA        NA
#> 2       P10         3     0.3303240    NA       NA       NA        NA
#> 3       P11         3     0.2469301    NA       NA       NA        NA
#> 4       P12         3     0.1029985    NA       NA       NA        NA
#> 5        P2         3     0.1738257    NA       NA       NA        NA
#> 6        P3         3     0.2105313    NA       NA       NA        NA
#> 7        P4         3     0.2653313    NA       NA       NA        NA
#> 8        P5         3     0.2447480    NA       NA       NA        NA
#> 9        P6         3     0.1725935    NA       NA       NA        NA
#> 10       P7         3     0.2267275    NA       NA       NA        NA
#> 11       P8         3     0.2954055    NA       NA       NA        NA
#> 12       P9         3     0.1714372    NA       NA       NA        NA
#>    alpha_raw_se alpha_raw_df alpha_covariates alpha_global alpha_global_se
#> 1            NA           NA                     0.8875006       0.1613071
#> 2            NA           NA                     0.8875006       0.1613071
#> 3            NA           NA                     0.8875006       0.1613071
#> 4            NA           NA                     0.8875006       0.1613071
#> 5            NA           NA                     0.8875006       0.1613071
#> 6            NA           NA                     0.8875006       0.1613071
#> 7            NA           NA                     0.8875006       0.1613071
#> 8            NA           NA                     0.8875006       0.1613071
#> 9            NA           NA                     0.8875006       0.1613071
#> 10           NA           NA                     0.8875006       0.1613071
#> 11           NA           NA                     0.8875006       0.1613071
#> 12           NA           NA                     0.8875006       0.1613071
#>     alpha_global_covariates alpha_prior_var alpha_prior_sd alpha_shrunk
#> 1  promoter;batch;replicate      0.02601998      0.1613071    0.8875006
#> 2  promoter;batch;replicate      0.02601998      0.1613071    0.8875006
#> 3  promoter;batch;replicate      0.02601998      0.1613071    0.8875006
#> 4  promoter;batch;replicate      0.02601998      0.1613071    0.8875006
#> 5  promoter;batch;replicate      0.02601998      0.1613071    0.8875006
#> 6  promoter;batch;replicate      0.02601998      0.1613071    0.8875006
#> 7  promoter;batch;replicate      0.02601998      0.1613071    0.8875006
#> 8  promoter;batch;replicate      0.02601998      0.1613071    0.8875006
#> 9  promoter;batch;replicate      0.02601998      0.1613071    0.8875006
#> 10 promoter;batch;replicate      0.02601998      0.1613071    0.8875006
#> 11 promoter;batch;replicate      0.02601998      0.1613071    0.8875006
#> 12 promoter;batch;replicate      0.02601998      0.1613071    0.8875006
#>    alpha_shrunk_se alpha_fixed_one alpha_diff_from_one
#> 1        0.1613071               1          -0.1124994
#> 2        0.1613071               1          -0.1124994
#> 3        0.1613071               1          -0.1124994
#> 4        0.1613071               1          -0.1124994
#> 5        0.1613071               1          -0.1124994
#> 6        0.1613071               1          -0.1124994
#> 7        0.1613071               1          -0.1124994
#> 8        0.1613071               1          -0.1124994
#> 9        0.1613071               1          -0.1124994
#> 10       0.1613071               1          -0.1124994
#> 11       0.1613071               1          -0.1124994
#> 12       0.1613071               1          -0.1124994
```

To reproduce the older fixed-ratio style response, `growth_exponent = 1`
gives the familiar `log2(luminescence / growth)` scale.

``` r

assay_fixed <- prepare_assay(
  screen,
  promoter = "promoter",
  compound = "compound",
  control = "DMSO",
  lux = "LUX.AUC_16",
  growth = "od_16h.measured",
  growth_exponent = 1,
  batch = "batch",
  replicate = "replicate"
)
```

## Fit the model workflow

[`fit_destress()`](https://muellsen.github.io/DStressR/reference/fit_destress.md)
fits promoter and compound effects while accounting for technical
covariates. The result table reports both the DMSO-relative total effect
and the promoter-specific effect after subtracting the compound-wide
effect.

``` r

fit <- fit_destress(
  assay,
  technical = c("batch", "replicate"),
  empirical_bayes = TRUE
)

tab <- results(fit)
tab <- adjust_pvalues(tab)

head(tab)
#>    promoter compound total_effect  total_se total_statistic total_pvalue
#> 1        P1       C1  -0.42759639 0.1264748      -3.3808826  0.001119764
#> 21       P1      C10   0.11008012 0.1264748       0.8703721  0.386701308
#> 34       P1      C11  -0.40523361 0.1264748      -3.2040665  0.001946616
#> 47       P1      C12   0.28971870 0.1264748       2.2907230  0.024611838
#> 58       P1      C13   0.04093903 0.1264748       0.3236932  0.747014983
#> 67       P1      C14  -0.06382475 0.1264748      -0.5046441  0.615196077
#>    additive_total_effect additive_total_se empty_vector_effect
#> 1                     NA                NA                  NA
#> 21                    NA                NA                  NA
#> 34                    NA                NA                  NA
#> 47                    NA                NA                  NA
#> 58                    NA                NA                  NA
#> 67                    NA                NA                  NA
#>    background_adjusted_effect global_effect  global_se global_statistic
#> 1                 -0.42759639   -0.54331614 0.03710069      -14.6443689
#> 21                 0.11008012    0.03412472 0.03710069        0.9197866
#> 34                -0.40523361   -0.53613209 0.03710069      -14.4507323
#> 47                 0.28971870    0.20219166 0.03710069        5.4498091
#> 58                 0.04093903    0.01268673 0.03710069        0.3419539
#> 67                -0.06382475   -0.01999624 0.03710069       -0.5389722
#>    global_pvalue low_rank_effect specific_effect specific_se specific_statistic
#> 1   2.391204e-24               0      0.11571976   0.1214185          0.9530654
#> 21  3.604494e-01               0      0.07595540   0.1214185          0.6255670
#> 34  5.196926e-24               0      0.13089849   0.1214185          1.0780771
#> 47  5.429453e-07               0      0.08752703   0.1214185          0.7208708
#> 58  7.332822e-01               0      0.02825231   0.1214185          0.2326854
#> 67  5.914028e-01               0     -0.04382851   0.1214185         -0.3609707
#>    specific_pvalue total_padj_global total_padj_by_promoter
#> 1        0.3434278        0.01414438             0.03016706
#> 21       0.5333795        0.61058101             0.53338111
#> 34       0.2842393        0.02123581             0.03016706
#> 47       0.4730904        0.11359310             0.09844735
#> 58       0.8166000        0.86820143             0.78633156
#> 67       0.7190729        0.79426626             0.70308123
#>    specific_padj_global specific_padj_by_promoter total_padj specific_padj
#> 1             0.9680909                 0.8211470 0.01414438     0.9680909
#> 21            0.9680909                 0.8211470 0.61058101     0.9680909
#> 34            0.9680909                 0.8211470 0.02123581     0.9680909
#> 47            0.9680909                 0.8211470 0.11359310     0.9680909
#> 58            0.9680909                 0.9332571 0.86820143     0.9680909
#> 67            0.9680909                 0.8716035 0.79426626     0.9680909
```

The key columns are:

- `total_effect`: DMSO-relative response for the promoter-compound pair.
- `global_effect`: average compound-wide response across promoters.
- `specific_effect`: promoter-specific deviation from the compound-wide
  effect.
- `specific_pvalue` and `specific_padj`: test and BH-adjusted p-value.

The direct fitting functions remain available for existing scripts, but
[`fit_destress()`](https://muellsen.github.io/DStressR/reference/fit_destress.md)
is the recommended entry point for new model-based analyses.

## Optional background reporter calibration

If a screen contains a matched background reporter, such as an Empty
Vector Control, the background reporter can be used during response
construction. When `background_promoter` is supplied, DStressR uses
Huber calibration by default; least-squares calibration and direct
subtraction remain available through `background_method`.

``` r

assay_bg <- prepare_assay(
  screen_with_evc,
  promoter = "promoter",
  compound = "compound",
  control = "DMSO",
  lux = "LUX.AUC_16",
  growth = "od_16h.measured",
  batch = "batch",
  replicate = "replicate",
  background_promoter = "EVC",
  background_by = c("compound", "batch", "replicate")
)

fit_bg <- fit_destress(
  assay_bg,
  technical = c("batch", "replicate"),
  empirical_bayes = TRUE
)
```

## Call hits

``` r

hits <- call_hits(
  tab,
  fdr = 0.05,
  lfc = 0.5,
  effect = "specific_effect",
  padj = "specific_padj"
)

table(hits$hit)
#> 
#> Downregulated        Not DE   Upregulated 
#>             5           472             3

head(
  hits[order(hits$specific_padj, -abs(hits$specific_effect)), ],
  10
)
#>     promoter compound total_effect  total_se total_statistic total_pvalue
#> 218      P12      C26   -1.3710711 0.1316085      -10.417805 1.490909e-16
#> 197      P11      C24   -1.6735292 0.1318927      -12.688569 7.533079e-21
#> 244      P10      C28    1.4895866 0.1379347       10.799216 2.736036e-17
#> 177      P10      C22   -1.6257162 0.1379347      -11.786129 3.610377e-19
#> 324      P10      C33    0.9705429 0.1379347        7.036249 6.043562e-10
#> 387      P11      C39   -1.2130464 0.1318927       -9.197224 3.599814e-14
#> 131       P5      C19    1.0692457 0.1415881        7.551803 6.070983e-11
#> 410       P8      C40   -0.7506326 0.1233521       -6.085283 3.799547e-08
#> 98        P8      C17   -0.2232267 0.1233521       -1.809670 7.410272e-02
#> 248       P4      C28   -0.2227294 0.1293115       -1.722425 8.885667e-02
#>     additive_total_effect additive_total_se empty_vector_effect
#> 218                    NA                NA                  NA
#> 197                    NA                NA                  NA
#> 244                    NA                NA                  NA
#> 177                    NA                NA                  NA
#> 324                    NA                NA                  NA
#> 387                    NA                NA                  NA
#> 131                    NA                NA                  NA
#> 410                    NA                NA                  NA
#> 98                     NA                NA                  NA
#> 248                    NA                NA                  NA
#>     background_adjusted_effect global_effect  global_se global_statistic
#> 218                 -1.3710711    0.11664730 0.03710069        3.1440739
#> 197                 -1.6735292   -0.20714467 0.03710069       -5.5833109
#> 244                  1.4895866    0.12623694 0.03710069        3.4025498
#> 177                 -1.6257162   -0.44954782 0.03710069      -12.1169675
#> 324                  0.9705429   -0.07827537 0.03710069       -2.1098092
#> 387                 -1.2130464   -0.23137805 0.03710069       -6.2364898
#> 131                  1.0692457    0.05542205 0.03710069        1.4938282
#> 410                 -0.7506326    0.01355452 0.03710069        0.3653441
#> 98                  -0.2232267    0.13065534 0.03710069        3.5216421
#> 248                 -0.2227294    0.12623694 0.03710069        3.4025498
#>     global_pvalue low_rank_effect specific_effect specific_se
#> 218  2.338259e-03               0      -1.4877184   0.1256741
#> 197  3.133514e-07               0      -1.4663846   0.1259103
#> 244  1.045062e-03               0       1.3633497   0.1309477
#> 177  8.648114e-20               0      -1.1761684   0.1309477
#> 324  3.800049e-02               0       1.0488183   0.1309477
#> 387  1.988889e-08               0      -0.9816684   0.1259103
#> 131  1.391544e-01               0       1.0138237   0.1340067
#> 410  7.158188e-01               0      -0.7641871   0.1188415
#> 98   7.115785e-04               0      -0.3538820   0.1188415
#> 248  1.045062e-03               0      -0.3489663   0.1237672
#>     specific_statistic specific_pvalue total_padj_global total_padj_by_promoter
#> 218         -11.837912    2.884547e-19      1.789091e-14           5.963636e-15
#> 197         -11.646263    6.628779e-19      3.615878e-18           3.013232e-19
#> 244          10.411407    1.534064e-16      4.377657e-15           5.472071e-16
#> 177          -8.981971    9.530509e-14      8.664905e-17           1.444151e-17
#> 324           8.009444    7.750641e-12      4.144157e-08           8.058082e-09
#> 387          -7.796569    2.022359e-11      3.455821e-12           7.199628e-13
#> 131           7.565471    5.710191e-11      4.856787e-09           2.428393e-09
#> 410          -6.430305    8.615587e-09      2.279728e-06           1.039530e-06
#> 98           -2.977764    3.841451e-03      2.355583e-01           2.571968e-01
#> 248          -2.819537    6.060215e-03      2.665700e-01           2.369511e-01
#>     specific_padj_global specific_padj_by_promoter   total_padj specific_padj
#> 218         1.384583e-16              1.153819e-17 1.789091e-14  1.384583e-16
#> 197         1.590907e-16              2.651511e-17 3.615878e-18  1.590907e-16
#> 244         2.454502e-14              6.136256e-15 4.377657e-15  2.454502e-14
#> 177         1.143661e-11              1.906102e-12 8.664905e-17  1.143661e-11
#> 324         7.440615e-10              1.033419e-10 4.144157e-08  7.440615e-10
#> 387         1.617887e-09              4.044719e-10 3.455821e-12  1.617887e-09
#> 131         3.915560e-09              2.284077e-09 4.856787e-09  3.915560e-09
#> 410         5.169352e-07              3.446235e-07 2.279728e-06  5.169352e-07
#> 98          2.048774e-01              7.682903e-02 2.355583e-01  2.048774e-01
#> 248         2.908903e-01              2.424086e-01 2.665700e-01  2.908903e-01
#>               hit
#> 218 Downregulated
#> 197 Downregulated
#> 244   Upregulated
#> 177 Downregulated
#> 324   Upregulated
#> 387 Downregulated
#> 131   Upregulated
#> 410 Downregulated
#> 98         Not DE
#> 248        Not DE
```

In real screens, the FDR cutoff should be paired with domain knowledge
and diagnostics. For DMSO-rich designs, empirical replicate/permutation
p-values can be used as an additional calibration check.

## Visualize the screen

The standard volcano plot labels the strongest promoter-compound pairs.

``` r

plot_volcano(
  hits,
  effect = "specific_effect",
  padj = "specific_padj",
  fdr = 0.05,
  lfc = 0.5,
  top_n = 8,
  top_promoters = 6
)
```

![](DStressR_files/figure-html/volcano-1.png)

A response heatmap gives a matrix-level view of the promoter-compound
response surface.

``` r

plot_response_heatmap(
  hits,
  value = "specific_effect",
  top_n_compounds = 30
)
```

![](DStressR_files/figure-html/heatmap-1.png)

For broader screens, clustered heatmaps help reveal compound and
promoter groups with similar response profiles.

``` r

plot_response_clustered_heatmap(
  hits,
  value = "specific_effect",
  top_n_compounds = 30,
  title = "Clustered DStressR response map"
)
```

![](DStressR_files/figure-html/clustered-1.png)

## Optional DGrowthR handoff

The current default hit model uses the exported growth summary column
supplied to
[`prepare_assay()`](https://muellsen.github.io/DStressR/reference/prepare_assay.md).
If growth curves have already been modeled with DGrowthR, the
DGrowthR-derived growth parameter can be joined explicitly before assay
preparation.

``` r

screen2 <- add_dgrowthr_growth(
  screen,
  object = dgrowthr_fit,
  by = "curve_id",
  model_covariate = "curve_id",
  growth_metric = "OD_16",
  output = "dgrowthr_od16"
)

assay2 <- prepare_assay(
  screen2,
  promoter = "promoter",
  compound = "compound",
  control = "DMSO",
  lux = "LUX.AUC_16",
  growth = "dgrowthr_od16",
  batch = "batch",
  replicate = "replicate"
)
```

This explicit handoff keeps the selected DStressR workflow reproducible,
while making DGrowthR-based sensitivity analyses straightforward.
