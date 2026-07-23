## R CMD check results

Local check:

```text
R CMD check DStressR_0.0.1.tar.gz --no-manual
Status: 1 NOTE
```

The local NOTE is due to the optional suggested package `nortest` not being
installed in the local check library. `nortest` is used only for the optional
Lilliefors normality test path.

The full manual check requires a local LaTeX installation providing
`pdflatex`.

## Package Scope

The CRAN package is intentionally lightweight. It includes package code,
documentation, tests, vignettes, and the public Binsfeld et al. reporter-screen
AUC/score data. Repository-only analysis workflows under `analysis/`,
`data-raw/`, and `paper/` are excluded from the source package with
`.Rbuildignore`.
