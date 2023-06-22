Simulation Results
================
Victor Tsang and David Warton
22 June, 2023

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

- MINMI point estimates aren’t as accurate as MLE under Uniform-Normal
  (UN) model
- MINMI point estimates appear to be more biased than MLE
- main issue is that MINMI has much higher sample variance - and not
  decreasing with sample size!!
- Coverage probability of MINMI is better than MLE for low sample size -
  but intervals very wide!
- The issue here is that sample minimum is not a good statistic for
  measurement error scenarios :(

#### Which value of `error_factor` do we want plots for?

``` r
error_fac_to_plot = 1 #change this value to look at plots for a different error_factor
```

------------------------------------------------------------------------

#### Load in the results

``` r
library(knitr)
library(tidyverse)
library(scales)
library(ggrepel)
library(gridExtra)
library(latex2exp)


load("data/synthetic-data.RData")
attach(synthetic.data.config)

RESULTS_PATH <- 'data/simResults-100-20230616.RData'
load(RESULTS_PATH)

head(results)
```

    ##   error_factor method    lower    point    upper point_runtime conf_int_runtime
    ## 1          0.0  MINMI 14830.29 14982.86 15016.26  7.891655e-05     7.891655e-05
    ## 2          0.0   UNci 14920.89 15017.52 15017.52            NA     5.306005e-02
    ## 3          0.0 UNwald 15017.52 15017.52 15017.52  0.000000e+00     7.192850e-03
    ## 4          0.5  MINMI 14833.54 15074.69 15630.77  1.979718e-01     1.979718e-01
    ## 5          0.5   UNci 14861.38 15026.43 15134.66            NA     1.551909e-01
    ## 6          0.5 UNwald 14899.04 15026.44 15153.83  6.499868e+01     6.108403e-02
    ##   B.lower B.point B.upper
    ## 1      NA      NA      NA
    ## 2     100     100     100
    ## 3     100     100     100
    ## 4     100     100     100
    ## 5     100     100     100
    ## 6     100     100     100

``` r
all_results=results
all_results$n.samples = synthetic.data.config$n.samples
all_results=all_results[0,]

n.samples=c(12,25,50,100)
for (iSample in 1:length(n.samples))
{
  RESULTS_PATH <- paste0("data/simResults-",n.samples[iSample],"-20230616.RData")
  load(RESULTS_PATH)
  all_results=tibble::add_row(
      all_results,
      error_factor = results$error_factor,
      method=results$method,
      lower=results$lower,
      point=results$point,
      upper=results$upper,
      point_runtime=results$point_runtime,
      conf_int_runtime=results$conf_int_runtime,
      B.lower=results$B.lower,
      B.point=results$B.point,
      B.upper=results$B.upper,
      n.samples=n.samples[iSample]
  )
}
```

``` r
all_results %>%
  group_by(method, error_factor, n.samples) %>%
  summarise(point.pct_na = mean(point,na.rm=TRUE),
            lower.pct_na = mean(lower,na.rm=TRUE),
            upper.pct_na = mean(upper,na.rm=TRUE))
```

    ## `summarise()` has grouped output by 'method', 'error_factor'. You can override
    ## using the `.groups` argument.

    ## # A tibble: 60 × 6
    ## # Groups:   method, error_factor [15]
    ##    method error_factor n.samples point.pct_na lower.pct_na upper.pct_na
    ##    <chr>         <dbl>     <dbl>        <dbl>        <dbl>        <dbl>
    ##  1 MINMI           0          12       15116.       13731.       15380.
    ##  2 MINMI           0          25       15062.       14434.       15192.
    ##  3 MINMI           0          50       15028.       14722.       15094.
    ##  4 MINMI           0         100       15014.       14863.       15048.
    ##  5 MINMI           0.5        12       15097.       13704.       15527.
    ##  6 MINMI           0.5        25       15021.       14364.       15421.
    ##  7 MINMI           0.5        50       14986.       14623.       15477.
    ##  8 MINMI           0.5       100       14957.       14714.       15515.
    ##  9 MINMI           1          12       15051.       13615.       15740.
    ## 10 MINMI           1          25       14968.       14228.       15758.
    ## # … with 50 more rows

# Point Estimates

#### Calculate Metrics

``` r
performance.point <- all_results %>%
  filter(!is.na(point)) %>% filter(error_factor == error_fac_to_plot) %>%
  group_by(method, n.samples) %>%
  summarise(MSE_000 = mean((point - theta.true)^2,na.rm=TRUE)/1000,
            bias = mean(point,na.rm=TRUE)-theta.true,
            variance_000 = var(point,na.rm=TRUE)/1000,
            avg_runtime = round(mean(point_runtime,na.rm=TRUE), 5))
```

    ## `summarise()` has grouped output by 'method'. You can override using the
    ## `.groups` argument.

``` r
performance.point.tbl = vector(mode = "list", length(n.samples))

for (i in 1:length(n.samples)) {
  performance.point.tbl[[i]] <- performance.point %>%
    filter(n.samples == n.samples[i]) %>%
    ungroup() %>%
    mutate(across(!c(method, n.samples, avg_runtime), round)) %>%
    mutate(avg_runtime = round(avg_runtime, digits = 5)) %>%
    arrange(MSE_000)
}

performance.point.tbl[[1]]
```

    ## # A tibble: 3 × 6
    ##   method n.samples MSE_000  bias variance_000 avg_runtime
    ##   <chr>      <dbl>   <dbl> <dbl>        <dbl>       <dbl>
    ## 1 MINMI         12     228    51          226      0.0801
    ## 2 UNci          12     240   205          199    NaN     
    ## 3 UNwald        12     240   205          199    303.

``` r
performance.point.tbl[[2]]
```

    ## # A tibble: 3 × 6
    ##   method n.samples MSE_000  bias variance_000 avg_runtime
    ##   <chr>      <dbl>   <dbl> <dbl>        <dbl>       <dbl>
    ## 1 UNci          25      77    81           70     NaN    
    ## 2 UNwald        25      77    81           70     228.   
    ## 3 MINMI         25     107   -32          106       0.146

``` r
performance.point.tbl[[3]]
```

    ## # A tibble: 3 × 6
    ##   method n.samples MSE_000  bias variance_000 avg_runtime
    ##   <chr>      <dbl>   <dbl> <dbl>        <dbl>       <dbl>
    ## 1 UNci          50      32    48           30     NaN    
    ## 2 UNwald        50      32    48           30     158.   
    ## 3 MINMI         50     184   -77          179       0.181

``` r
performance.point.tbl[[4]]
```

    ## # A tibble: 3 × 6
    ##   method n.samples MSE_000  bias variance_000 avg_runtime
    ##   <chr>      <dbl>   <dbl> <dbl>        <dbl>       <dbl>
    ## 1 UNci         100      15    32           14     NaN    
    ## 2 UNwald       100      15    32           14     106.   
    ## 3 MINMI        100     198  -104          187       0.559

Ignore run-times - I used that to store SEs for UN

#### Pivot to make plots and filter to `error_factor`

``` r
performance.point.long <- performance.point %>%
  rename(Method = method, n.samples = n.samples, Bias = bias, Var_000 = variance_000, Runtime = avg_runtime) %>%
  pivot_longer(cols=c(MSE_000, Bias, Var_000, Runtime), names_to = "Metric")
  
performance.point.long
```

    ## # A tibble: 48 × 4
    ## # Groups:   Method [3]
    ##    Method n.samples Metric     value
    ##    <chr>      <dbl> <chr>      <dbl>
    ##  1 MINMI         12 MSE_000 228.    
    ##  2 MINMI         12 Bias     50.9   
    ##  3 MINMI         12 Var_000 226.    
    ##  4 MINMI         12 Runtime   0.0801
    ##  5 MINMI         25 MSE_000 107.    
    ##  6 MINMI         25 Bias    -31.9   
    ##  7 MINMI         25 Var_000 106.    
    ##  8 MINMI         25 Runtime   0.146 
    ##  9 MINMI         50 MSE_000 184.    
    ## 10 MINMI         50 Bias    -76.6   
    ## # … with 38 more rows

### Plots

``` r
metrics = unique(performance.point.long$Metric)
performance.point_estimates.plots = lapply(metrics,
  function(met) {
    p = ggplot(data = filter(performance.point.long, Metric == met),
               mapping = aes(x = n.samples, y = value, colour = reorder(Method, value, decreasing=T))) +
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
                                    "UNci" = "darkblue"))
    
    if (met %in% c("MSE", "Runtime")) {
      p = p+scale_y_log10(labels = label_comma())
    }
    p
  }
)

performance.point_estimates.plots[[1]] = performance.point_estimates.plots[[1]] + ylab("000's")
performance.point_estimates.plots[[2]] = performance.point_estimates.plots[[2]] + ylab("Years")
performance.point_estimates.plots[[3]] = performance.point_estimates.plots[[3]] + ylab("000's")

performance.point_estimates.plots[[1]]
```

![](ResultsForNSimulations_files/figure-gfm/unnamed-chunk-5-1.png)<!-- -->

``` r
performance.point_estimates.plots[[2]]
```

![](ResultsForNSimulations_files/figure-gfm/unnamed-chunk-5-2.png)<!-- -->

``` r
performance.point_estimates.plots[[3]]
```

![](ResultsForNSimulations_files/figure-gfm/unnamed-chunk-5-3.png)<!-- -->

## Commentary

1.  MSE:
    1.  MINMI has much higher MSE than UN-MLE
    2.  MINMI MSE does not decrease as sample size increases (for large
        n)
2.  Bias:
    1.  MINMI is more biased than UN (mean vs median thing I guess)
3.  Variance:
    1.  MINMI estimates generally have more variance than UN, and not
        decreasing with sample size (for large n)

``` r
performance.point_estimates.plot.grid = do.call(grid.arrange, performance.point_estimates.plots)
```

    ## Warning: Removed 4 row(s) containing missing values (geom_path).

    ## Warning: Removed 4 rows containing missing values (geom_point).

![](ResultsForNSimulations_files/figure-gfm/unnamed-chunk-6-1.png)<!-- -->

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
performance.CI <- all_results %>%
  filter(!is.na(conf_int_runtime)) %>% filter(error_factor == error_fac_to_plot) %>%
  mutate(width = upper - lower,
         contains_theta = ifelse(theta.true > lower & theta.true < upper, 1, 0)) %>%
  group_by(n.samples, method) %>%
  summarise(Coverage = round(mean(contains_theta, na.rm=TRUE) * 100, 1),
            `Average Width` = round(mean(width, na.rm=TRUE), 2),
            `Average Runtime` = round(mean(conf_int_runtime, na.rm=TRUE), 5)) %>%
  ungroup() %>%
  arrange(method, n.samples)
```

    ## `summarise()` has grouped output by 'n.samples'. You can override using the
    ## `.groups` argument.

``` r
performance.CI.long <- performance.CI %>%
  rename(n.samples = n.samples, Method = method, Width = `Average Width`, Runtime = `Average Runtime`) %>%
  pivot_longer(cols=c(Coverage, Width, Runtime),
               names_to = "Metric")
  
performance.CI.long
```

    ## # A tibble: 36 × 4
    ##    n.samples Method Metric       value
    ##        <dbl> <chr>  <chr>        <dbl>
    ##  1        12 MINMI  Coverage   93.9   
    ##  2        12 MINMI  Width    2125.    
    ##  3        12 MINMI  Runtime     0.0801
    ##  4        25 MINMI  Coverage   96.7   
    ##  5        25 MINMI  Width    1530.    
    ##  6        25 MINMI  Runtime     0.146 
    ##  7        50 MINMI  Coverage   94.5   
    ##  8        50 MINMI  Width    1566.    
    ##  9        50 MINMI  Runtime     0.181 
    ## 10       100 MINMI  Coverage   94     
    ## # … with 26 more rows

## Coverage Probability

``` r
conf_int.coverage.plot <- performance.CI.long %>%
  filter(Metric == "Coverage") %>%
  ggplot(aes(x=n.samples, y=value, colour=reorder(Method, value, decreasing=T))) +
  geom_point() +
  geom_line(linewidth=0.5) +
  geom_label_repel(aes(label = value)) +
  theme_bw() +
  labs(y = "Coverage probability", colour="Method", title="Coverage Probabilities") +
  scale_y_continuous(breaks=c(0, 25, 50, 75, 95, 100)) +
  theme(rect = element_rect(fill = "transparent")) +
  scale_color_manual(values = c("GRIWM" = "#F8766D", "GRIWM-corrected" = "#619CFF", "MINMI" = "#00BA38", "UNci" = "darkblue", "UNwald"="red"))
```

    ## Warning: Ignoring unknown parameters: linewidth

``` r
conf_int.coverage.plot
```

![](ResultsForNSimulations_files/figure-gfm/unnamed-chunk-8-1.png)<!-- -->

## Widths

``` r
conf_int.width.plot <- performance.CI.long %>%
  filter(Metric == "Width") %>%
  ggplot(aes(x=n.samples, y=value, colour=reorder(Method, value, decreasing=T))) +
  geom_point() +
  geom_line(linewidth=0.5) +
  theme_bw() +
  labs(y = "Years", colour="Method", title="Average Width of Estimated Confidence Intervals") +
  theme(rect = element_rect(fill = "transparent")) +
  scale_color_manual(values = c("GRIWM" = "#F8766D", "GRIWM-corrected" = "#619CFF", "MINMI" = "#00BA38", "UNci"="darkblue", "UNwald"="red"))
```

    ## Warning: Ignoring unknown parameters: linewidth

``` r
conf_int.width.plot
```

![](ResultsForNSimulations_files/figure-gfm/unnamed-chunk-9-1.png)<!-- -->

## Runtime

``` r
conf_int.runtime.plot <- performance.CI.long %>%
  filter(Metric == "Runtime") %>%
  ggplot(aes(x=n.samples, y=value, colour=reorder(Method, value, decreasing=T))) +
  geom_point() +
  geom_line(linewidth=0.5) +
  theme_bw() +
  scale_y_continuous(trans=scales::log10_trans()) +
  labs(y = "Seconds", colour="Method", title="Average Runtime of Confidence Interval Estimation") +
  theme(rect = element_rect(fill = "transparent")) +
  scale_color_manual(values = c("GRIWM" = "#F8766D", "GRIWM-corrected" = "#619CFF", "MINMI" = "#00BA38", "UNci"="darkblue", "UNwald"="red"))
```

    ## Warning: Ignoring unknown parameters: linewidth

``` r
conf_int.runtime.plot
```

![](ResultsForNSimulations_files/figure-gfm/unnamed-chunk-10-1.png)<!-- -->

## Commentary

1.  Coverage Probability:
    1.  MINMI generally has good coverage probability
    2.  UN-MLE methods have good coverage at moderate-large sample sizes
        (but not for small n)
    3.  UN-Wald is slower to converge to desired coverage probability
        (symmetric CI)
2.  Confidence Interval Widths:
    1.  MINMI has way wider CIs
    2.  CI width does not decrease as sample size increases
3.  Runtime
    1.  MINMI takes longer! Although I was using B=100, might have been
        fairer to compare at B=100 instead of doing choose B first

#### Extra bonus: is the sampling distribution of MLE Gaussian for ![\sigma\>0](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Csigma%3E0 "\sigma>0")?

``` r
errors=unique(results$error_factor)
nError=length(errors)
par(mfrow=c(5,3),mgp=c(1.75,0.75,0),mar=c(3,2,0,0),oma=c(0,2,2,0))
for(iError in 1:nError)
{
  tmp=results%>%filter(method=="UNci" & error_factor==errors[iError]) %>% select(lower,point,upper)
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

![](ResultsForNSimulations_files/figure-gfm/distTheta-1.png)<!-- -->

Um, yes! Not at , as expected, because this is a sample minimum.

Quantiles also seem to be approx normal with no outliers (except at
where it looks like there is some non-convergence).

#### How good are standard error estimates?

``` r
errors=unique(results$error_factor)
nError=length(errors)
par(mfrow=c(3,2),mgp=c(1.75,0.75,0),mar=c(3,2,0,0),oma=c(0,2,2,0))
for(iError in 1:nError)
{
  tmp=results%>%filter(method=="UNwald" & error_factor==errors[iError]) %>% select(point,point_runtime)
  nError=length(errors)
  hist(tmp$point_runtime,xlab="SE",ylab="",main="")
  mtext(paste0("error_fac=",errors[iError]),2,line=2,font=2,cex=0.8)
  abline(v=sd(tmp$point),col="red")
}
```

![](ResultsForNSimulations_files/figure-gfm/waldSE-1.png)<!-- -->
