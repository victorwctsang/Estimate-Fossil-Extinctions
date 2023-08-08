Simulation Results
================
Victor Tsang
08 August, 2023

- <a href="#tldr" id="toc-tldr">TL;DR</a>
- <a href="#point-estimates" id="toc-point-estimates">Point Estimates</a>
  - <a href="#plots" id="toc-plots">Plots</a>
  - <a href="#commentary" id="toc-commentary">Commentary</a>
- <a href="#confidence-intervals" id="toc-confidence-intervals">Confidence
  Intervals</a>
  - <a href="#coverage-probability" id="toc-coverage-probability">Coverage
    Probability</a>
  - <a href="#widths" id="toc-widths">Widths</a>
  - <a href="#runtime" id="toc-runtime">Runtime</a>
  - <a href="#commentary-1" id="toc-commentary-1">Commentary</a>

# TL;DR

- MINMI point estimates aren’t as accurate as other methods (MLE) in
  high measurement error variation scenarios, because the minimum is not
  sufficient in the measurement error setting
- MINMI point estimates appear to be more biased in high measurement
  error settings and also more variable.
- As expected, MLE_INV is much slower than MINMI and asymptotic MLE
  approaches.
- Asymptotic MLE approaches have poor coverage probability when sample
  size and measurement error are both small, and as expected, MINMI and
  MLE_INV do fine on coverage.
- MINMI has wide confidence intervals when there is measurement error
  (because of inefficiency using minimum as statistic)

------------------------------------------------------------------------

#### Load in the results

``` r
library(knitr)
library(tidyverse)
library(scales)
library(ggrepel)
library(gridExtra)
library(latex2exp)


load("data/synthetic-data-24-20230808.RData")
attach(synthetic.data.config)

RESULTS_PATH <- 'data/simResults-24-20230808.RData'
load(RESULTS_PATH)

head(results)
```

    ##   which_sim n.samples error_factor       method     lower    point    upper
    ## 1         1        24          0.0        MINMI  8854.023 10161.96 10431.94
    ## 2         1        24          0.0    mlereginv  8798.163       NA 10429.64
    ## 3         1        24          0.0   reginvUNci  9645.661 10442.03 10442.03
    ## 4         1        24          0.0 reginvUNwald 10442.030 10442.03 10442.03
    ## 5         2        24          0.5        MINMI  8849.419 10158.06 10461.70
    ## 6         2        24          0.5    mlereginv  8687.795       NA 10446.43
    ##   point_runtime conf_int_runtime B.lower B.point B.upper
    ## 1  0.0003247261     3.247261e-04      NA      NA      NA
    ## 2  0.0001287460     1.849224e+01      NA      NA      NA
    ## 3  0.0000000000     1.747608e-03      NA      NA      NA
    ## 4  0.0000000000     1.733303e-04      NA      NA      NA
    ## 5  0.0539743900     5.397439e-02     100     100     100
    ## 6  0.0010674000     2.136510e+01      NA      NA      NA

``` r
results %>%
  group_by(method, error_factor) %>%
  summarise(point.pct_na = mean(point,na.rm=TRUE),
            lower.pct_na = mean(lower,na.rm=TRUE),
            upper.pct_na = mean(upper,na.rm=TRUE))
```

    ## `summarise()` has grouped output by 'method'. You can override using the
    ## `.groups` argument.

    ## # A tibble: 24 × 5
    ## # Groups:   method [4]
    ##    method    error_factor point.pct_na lower.pct_na upper.pct_na
    ##    <chr>            <dbl>        <dbl>        <dbl>        <dbl>
    ##  1 MINMI              0         10134.        8822.       10405.
    ##  2 MINMI              0.5       10117.        8804.       10423.
    ##  3 MINMI              1         10129.        8817.       10501.
    ##  4 MINMI              2         10069.        8740.       10595.
    ##  5 MINMI              4         10049.        8639.       10888.
    ##  6 MINMI              8          9908.        8173.       11345.
    ##  7 mlereginv          0           NaN         8825.       10404.
    ##  8 mlereginv          0.5         NaN         8803.       10423.
    ##  9 mlereginv          1           NaN         8822.       10499.
    ## 10 mlereginv          2           NaN         8742.       10583.
    ## # … with 14 more rows

# Point Estimates

#### Calculate Metrics

``` r
performance.point <- results %>%
  filter(!is.na(point)) %>%
  group_by(error_factor, method) %>%
  summarise(MSE_000 = mean((point - theta.true)^2,na.rm=TRUE)/1000,
            bias = mean(point,na.rm=TRUE)-theta.true,
            variance_000 = var(point,na.rm=TRUE)/1000,
            avg_runtime = round(mean(point_runtime,na.rm=TRUE), 5))
```

    ## `summarise()` has grouped output by 'error_factor'. You can override using the
    ## `.groups` argument.

``` r
performance.point.tbl = vector(mode = "list", length(error_factors))

for (i in 1:length(error_factors)) {
  performance.point.tbl[[i]] <- performance.point %>%
    filter(error_factor == error_factors[i]) %>%
    ungroup() %>%
    mutate(across(!c(error_factor, method, avg_runtime), round)) %>%
    mutate(avg_runtime = round(avg_runtime, digits = 5)) %>%
    arrange(MSE_000)
}

performance.point.tbl[[1]]
```

    ## # A tibble: 3 × 6
    ##   error_factor method       MSE_000  bias variance_000 avg_runtime
    ##          <dbl> <chr>          <dbl> <dbl>        <dbl>       <dbl>
    ## 1            0 MINMI            198   134          181      0.0001
    ## 2            0 reginvUNci       342   415          170      0     
    ## 3            0 reginvUNwald     342   415          170      0

``` r
performance.point.tbl[[2]]
```

    ## # A tibble: 3 × 6
    ##   error_factor method       MSE_000  bias variance_000 avg_runtime
    ##          <dbl> <chr>          <dbl> <dbl>        <dbl>       <dbl>
    ## 1          0.5 MINMI            181   117          167      0.0445
    ## 2          0.5 reginvUNci       261   319          159    108.    
    ## 3          0.5 reginvUNwald     261   319          159    108.

``` r
performance.point.tbl[[3]]
```

    ## # A tibble: 3 × 6
    ##   error_factor method       MSE_000  bias variance_000 avg_runtime
    ##          <dbl> <chr>          <dbl> <dbl>        <dbl>       <dbl>
    ## 1            1 MINMI            178   129          162      0.0474
    ## 2            1 reginvUNci       241   294          155    169.    
    ## 3            1 reginvUNwald     241   294          155    169.

``` r
performance.point.tbl[[4]]
```

    ## # A tibble: 3 × 6
    ##   error_factor method       MSE_000  bias variance_000 avg_runtime
    ##          <dbl> <chr>          <dbl> <dbl>        <dbl>       <dbl>
    ## 1            2 MINMI            191    69          187      0.0541
    ## 2            2 reginvUNci       221   206          179    264.    
    ## 3            2 reginvUNwald     221   206          179    264.

``` r
performance.point.tbl[[5]]
```

    ## # A tibble: 3 × 6
    ##   error_factor method       MSE_000  bias variance_000 avg_runtime
    ##          <dbl> <chr>          <dbl> <dbl>        <dbl>       <dbl>
    ## 1            4 MINMI            293    49          291      0.0589
    ## 2            4 reginvUNci       298   181          265    410.    
    ## 3            4 reginvUNwald     298   181          265    410.

#### Pivot to make plots

``` r
performance.point.long <- performance.point %>%
  rename(Error = error_factor, Method = method, Bias = bias, Var_000 = variance_000, Runtime = avg_runtime) %>%
  pivot_longer(cols=c(MSE_000, Bias, Var_000, Runtime), names_to = "Metric")
  
performance.point.long
```

    ## # A tibble: 72 × 4
    ## # Groups:   Error [6]
    ##    Error Method       Metric     value
    ##    <dbl> <chr>        <chr>      <dbl>
    ##  1     0 MINMI        MSE_000 198.    
    ##  2     0 MINMI        Bias    134.    
    ##  3     0 MINMI        Var_000 181.    
    ##  4     0 MINMI        Runtime   0.0001
    ##  5     0 reginvUNci   MSE_000 342.    
    ##  6     0 reginvUNci   Bias    415.    
    ##  7     0 reginvUNci   Var_000 170.    
    ##  8     0 reginvUNci   Runtime   0     
    ##  9     0 reginvUNwald MSE_000 342.    
    ## 10     0 reginvUNwald Bias    415.    
    ## # … with 62 more rows

### Plots

``` r
metrics = unique(performance.point.long$Metric)
performance.point_estimates.plots = lapply(metrics,
  function(met) {
    p = ggplot(data = filter(performance.point.long, Metric == met),
               mapping = aes(x = Error, y = value, colour = reorder(Method, value, decreasing=T))) +
      geom_line() +
      geom_point() +
      theme_bw() +
      labs(title = paste(met, "by Error"), ylab=NULL, colour = "Method") +
      theme(rect = element_rect(fill = "transparent")) +
      scale_color_manual(values = c("MINMI" = "#00BA38",
                                    "MLE" = "#619CFF",
                                    "BA-MLE" = "purple",
                                    "Strauss" = "orange",
                                    "GRIWM-corrected" = "darkgray",
                                    "GRIWM" = "maroon",
                                    "reginvUNci" = "darkblue",
                                    "reginvUNwald" = "red",
                                    "mlereginv" = "purple",
                                    "mleInvAW" = "pink"))
    
    if (met %in% c("MSE", "Runtime")) {
      p = p+scale_y_log10(labels = label_comma())
    }
    p
  }
)

performance.point_estimates.plots[[1]] = performance.point_estimates.plots[[1]] + ylab("000's")
performance.point_estimates.plots[[2]] = performance.point_estimates.plots[[2]] + ylab("Years")
performance.point_estimates.plots[[3]] = performance.point_estimates.plots[[3]] + ylab("000's")
performance.point_estimates.plots[[4]] = performance.point_estimates.plots[[4]] + ylab("Seconds")

performance.point_estimates.plots[[1]]
```

![](ResultsSimulations_files/figure-gfm/unnamed-chunk-6-1.png)<!-- -->

``` r
performance.point_estimates.plots[[2]]
```

![](ResultsSimulations_files/figure-gfm/unnamed-chunk-6-2.png)<!-- -->

``` r
performance.point_estimates.plots[[3]]
```

![](ResultsSimulations_files/figure-gfm/unnamed-chunk-6-3.png)<!-- -->

``` r
performance.point_estimates.plots[[4]]
```

    ## Warning: Transformation introduced infinite values in continuous y-axis
    ## Transformation introduced infinite values in continuous y-axis

![](ResultsSimulations_files/figure-gfm/unnamed-chunk-6-4.png)<!-- -->

## Commentary

1.  MSE:
    1.  MINMI generally produces estimates with similar MSE to the MLE
    2.  MINMI had the worst MSE in
        ![4\sigma](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;4%5Csigma "4\sigma")
        scenarios and was moderately bad in the
        ![0\sigma](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;0%5Csigma "0\sigma")
        scenario
2.  Bias:
    1.  MINMI does OK but has poor bias in the
        ![4\sigma](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;4%5Csigma "4\sigma")
        scenario, when minimum is clearly not the optimal statistic.
3.  Variance:
    1.  MINMI estimates generally have more variance than the other
        methods, especially in high measurement error scenarios.
4.  Runtime:
    1.  MINMI is comparable to MLE

``` r
performance.point_estimates.plot.grid = do.call(grid.arrange, performance.point_estimates.plots)
```

    ## Warning: Transformation introduced infinite values in continuous y-axis
    ## Transformation introduced infinite values in continuous y-axis

![](ResultsSimulations_files/figure-gfm/unnamed-chunk-7-1.png)<!-- -->

``` r
performance.point_estimates.plot.grid
```

    ## TableGrob (2 x 2) "arrange": 4 grobs
    ##   z     cells    name           grob
    ## 1 1 (1-1,1-1) arrange gtable[layout]
    ## 2 2 (1-1,2-2) arrange gtable[layout]
    ## 3 3 (2-2,1-1) arrange gtable[layout]
    ## 4 4 (2-2,2-2) arrange gtable[layout]

# Confidence Intervals

#### Calculate Metrics and Pivot

``` r
performance.CI <- results %>%
  filter(!is.na(conf_int_runtime)) %>%
  mutate(width = upper - lower,
         contains_theta = ifelse(theta.true > lower & theta.true < upper, 1, 0)) %>%
#         contains_theta = ifelse(theta.true > lower, 1, 0)) %>%
#         contains_theta = ifelse(theta.true < upper, 1, 0)) %>%
  group_by(error_factor, method) %>%
  summarise(Coverage = round(mean(contains_theta, na.rm=TRUE) * 100, 1),
            `Average Width` = round(mean(width, na.rm=TRUE), 2),
            `Average Runtime` = round(mean(conf_int_runtime, na.rm=TRUE), 5)) %>%
  ungroup() %>%
  arrange(method, error_factor)
```

    ## `summarise()` has grouped output by 'error_factor'. You can override using the
    ## `.groups` argument.

``` r
performance.CI.long <- performance.CI %>%
  rename(Error = error_factor, Method = method, Width = `Average Width`, Runtime = `Average Runtime`) %>%
  pivot_longer(cols=c(Coverage, Width, Runtime),
               names_to = "Metric")
  
performance.CI.long
```

    ## # A tibble: 72 × 4
    ##    Error Method Metric       value
    ##    <dbl> <chr>  <chr>        <dbl>
    ##  1   0   MINMI  Coverage   94.5   
    ##  2   0   MINMI  Width    1582.    
    ##  3   0   MINMI  Runtime     0.0001
    ##  4   0.5 MINMI  Coverage   95.4   
    ##  5   0.5 MINMI  Width    1619.    
    ##  6   0.5 MINMI  Runtime     0.0445
    ##  7   1   MINMI  Coverage   95.6   
    ##  8   1   MINMI  Width    1683.    
    ##  9   1   MINMI  Runtime     0.0474
    ## 10   2   MINMI  Coverage   95     
    ## # … with 62 more rows

## Coverage Probability

``` r
conf_int.coverage.plot <- performance.CI.long %>%
  filter(Metric == "Coverage") %>%
  ggplot(aes(x=Error, y=value, colour=reorder(Method, value, decreasing=T))) +
  geom_point() +
  geom_line(linewidth=0.5) +
  geom_label_repel(aes(label = value)) +
  theme_bw() +
  labs(y = "Years", colour="Method", title="Coverage Probabilities") +
  scale_y_continuous(breaks=c(0, 25, 50, 75, 95, 100)) +
  theme(rect = element_rect(fill = "transparent")) +
  scale_color_manual(values = c("GRIWM" = "#F8766D", "GRIWM-corrected" = "#619CFF", "MINMI" = "#00BA38",                                     "reginvUNci" = "darkblue",
                                    "reginvUNwald" = "red",
                                    "mlereginv" = "purple",
                                    "mleInvAW" = "pink"))
```

    ## Warning: Ignoring unknown parameters: linewidth

``` r
conf_int.coverage.plot
```

    ## Warning: ggrepel: 7 unlabeled data points (too many overlaps). Consider
    ## increasing max.overlaps

![](ResultsSimulations_files/figure-gfm/unnamed-chunk-9-1.png)<!-- -->

## Widths

``` r
conf_int.width.plot <- performance.CI.long %>%
  filter(Metric == "Width") %>%
  ggplot(aes(x=Error, y=value, colour=reorder(Method, value, decreasing=T))) +
  geom_point() +
  geom_line(linewidth=0.5) +
  theme_bw() +
  labs(y = "Years", colour="Method", title="Average Width of Estimated Confidence Intervals") +
  theme(rect = element_rect(fill = "transparent")) +
  scale_color_manual(values = c("GRIWM" = "#F8766D", "GRIWM-corrected" = "#619CFF", "MINMI" = "#00BA38",
                                                                    "reginvUNci" = "darkblue",
                                    "reginvUNwald" = "red",
                                    "mlereginv" = "purple",
                                    "mleInvAW" = "pink"))
```

    ## Warning: Ignoring unknown parameters: linewidth

``` r
conf_int.width.plot
```

![](ResultsSimulations_files/figure-gfm/unnamed-chunk-10-1.png)<!-- -->

## Runtime

``` r
conf_int.runtime.plot <- performance.CI.long %>%
  filter(Metric == "Runtime") %>%
  ggplot(aes(x=Error, y=value, colour=reorder(Method, value, decreasing=T))) +
  geom_point() +
  geom_line(linewidth=0.5) +
  theme_bw() +
  scale_y_continuous(trans=scales::log10_trans()) +
  labs(y = "Seconds", colour="Method", title="Average Runtime of Confidence Interval Estimation") +
  theme(rect = element_rect(fill = "transparent")) +
  scale_color_manual(values = c("GRIWM" = "#F8766D", "GRIWM-corrected" = "#619CFF", "MINMI" = "#00BA38", 
                                                                    "reginvUNci" = "darkblue",
                                    "reginvUNwald" = "red",
                                    "mlereginv" = "purple",
                                    "mleInvAW" = "pink"))
```

    ## Warning: Ignoring unknown parameters: linewidth

``` r
conf_int.runtime.plot
```

![](ResultsSimulations_files/figure-gfm/unnamed-chunk-11-1.png)<!-- -->

## Commentary

1.  Coverage Probability:
    1.  MINMI and MLE_INV have good coverage as expected
    2.  Asymptotic MLE methods have poor coverage for small n and small
        measurement error, especially Wald
2.  Confidence Interval Widths:
    1.  MINMI has consistently wider CI’s
    2.  MLE_INV is also a bit wider, expected since other methods have
        undercoverage, although diff bigger than expected
3.  Runtime
    1.  MLE_INV much slower than everything else, and surprisingly,
        asymptotic MLE methods faster than MINMI.

#### Bonus: measurement error variation relative to our sampling error variation?

``` r
pct_sigma_sampling <- 4*fossil.sd / (K-theta.true)

tibble(index = 1:n.samples, pct_sigma_sampling) %>%
  mutate(label = ifelse(pct_sigma_sampling > 0.3, percent(pct_sigma_sampling), "")) %>%
  ggplot(aes(x=index, y=pct_sigma_sampling)) +
  geom_point() +
  geom_label_repel(aes(label=label)) + 
  labs(x = 'Sample Index', y = '% of K - theta', title="Measurement Error Variation Relative to Sample Error Variation", subtitle = "(4sigma scenario)") +
  scale_y_continuous(labels = percent_format())
```

![](ResultsSimulations_files/figure-gfm/unnamed-chunk-12-1.png)<!-- -->

``` r
tibble(index = 1:n.samples, pct_sigma_sampling) %>%
  ggplot(aes(x=pct_sigma_sampling)) +
  geom_histogram(binwidth=0.05) +
  scale_x_continuous(labels = percent_format())
```

![](ResultsSimulations_files/figure-gfm/unnamed-chunk-13-1.png)<!-- -->

Under
![4\sigma](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;4%5Csigma "4\sigma")
scenario, we have a right skewed distribution. Our fossils are mostly
\<30% of
![K-\theta](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;K-%5Ctheta "K-\theta"),
but we do get some samples with super large measurement error variation.
Perhaps these cause problems?

#### Extra bonus: is the sampling distribution of MLE Gaussian for ![\sigma\>0](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Csigma%3E0 "\sigma>0")?

``` r
errors=unique(results$error_factor)
nError=length(errors)
par(mfrow=c(5,3),mgp=c(1.75,0.75,0),mar=c(3,2,0,0),oma=c(0,2,2,0))
for(iError in 1:nError)
{
  tmp=results%>%filter(method=="reginvUNci" & error_factor==errors[iError]) %>% select(lower,point,upper)
  hist(tmp$lower,xlab="theta_lower",ylab="",main="")
  mtext(paste0("error_fac=",errors[iError]),2,line=2,font=2,cex=0.8)
  if(iError==1)
    mtext("theta_lower",3,font=2,cex=0.8)
  hist(tmp$point,xlab="theta_hat",ylab="",main="")
  if(iError==1)
    mtext("theta_hat",3,font=2)
  hist(tmp$point,xlab="theta_upper",ylab="",main="")
  if(iError==1)
    mtext("theta_upper",3,font=2)
}
```

![](ResultsSimulations_files/figure-gfm/distTheta-1.png)<!-- -->![](ResultsSimulations_files/figure-gfm/distTheta-2.png)<!-- -->

Um, yes! Not at , as expected, because this is a sample minimum.

Quantiles also seem to be approx normal with no outliers (except at 0=4
where it looks like there is some non-convergence).

#### How good are standard error estimates?

``` r
errors=unique(results$error_factor)
nError=length(errors)
par(mfrow=c(3,2),mgp=c(1.75,0.75,0),mar=c(3,2,0,0),oma=c(0,2,2,0))
for(iError in 1:nError)
{
  tmp=results%>%filter(method=="reginvUNwald" & error_factor==errors[iError]) %>% select(point,point_runtime)
  nError=length(errors)
  hist(tmp$point_runtime,xlab="SE",ylab="",main="")
  abline(v=sd(tmp$point),col="red")
}
```

![](ResultsSimulations_files/figure-gfm/waldSE-1.png)<!-- -->
