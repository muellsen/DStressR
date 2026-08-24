## R CMD check results

Local check:

```text
_R_CHECK_FORCE_SUGGESTS_=false R CMD check DStressR_0.0.1.tar.gz
Status: OK
```

The package builds and checks locally without errors, warnings, or notes. The
source package excludes repository-only analysis outputs, manuscript files, and
local VennDiagram log files via `.Rbuildignore`.

## Package Scope

The CRAN package is intentionally lightweight. It includes package code,
documentation, tests, vignettes, and the public Binsfeld et al. reporter-screen
AUC/score data. Repository-only analysis workflows under `analysis/`,
`data-raw/`, and `paper/` are excluded from the source package with
`.Rbuildignore`.
