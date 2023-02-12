Simulation Results
================
Victor Tsang
13 February, 2023

``` r
RESULTS_PATH <- 'data/simResults.RData'
load(RESULTS_PATH)
head(results)
```

    ##   error_factor          method    lower    point    upper point_runtime
    ## 1            0         Strauss       NA 14973.11       NA  2.193451e-05
    ## 2            0             MLE       NA 15073.05       NA  1.435041e-03
    ## 3            0          BA-MLE       NA 14974.51       NA  2.193451e-05
    ## 4            0           MINMI 14695.81 15004.28 15070.56  2.781296e-02
    ## 5            0           GRIWM 14779.00 14779.00 14779.00  1.442542e+00
    ## 6            0 GRIWM-corrected 15004.00 15004.00 15004.00  3.216715e+01
    ##   conf_int_runtime B.lower B.point B.upper
    ## 1               NA      NA      NA      NA
    ## 2               NA      NA      NA      NA
    ## 3               NA      NA      NA      NA
    ## 4       0.02781296      NA      NA      NA
    ## 5       1.44254208      NA      NA      NA
    ## 6      32.16714787      NA      NA      NA

### Point Estimates

### Confidence Intervals
