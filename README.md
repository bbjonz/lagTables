# lagTables
Computes lag sequential statistics for categorical event sequence data and
renders the results as html tables (and plots, if requested), following the
log-linear approach to lag sequential analysis described in Poole (2000).

- `trprobs()` — observed/expected frequency and transition-probability
  tables with standardized residuals, optionally by group, with
  standardized-residual bar plots
- `lagmodels()` — log-linear models testing lag order/dependency
- `shmodels()` — stationarity test (splits each sequence in half and
  tests the lag × time interaction)

# Installation
```r
remotes::install_github("bbjonz/lagTables")
```
