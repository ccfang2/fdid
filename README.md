# fdid <img src="man/figures/fdid_badge.png" align="right" alt="" width="155" />

[![PRs Welcome](https://img.shields.io/badge/PRs-welcome-brightgreen.svg?style=flat-square)](https://makeapullrequest.com)
![GitHub last commit](https://img.shields.io/github/last-commit/ccfang2/fdid?logo=GitHub)
![GitHub repo size](https://img.shields.io/github/repo-size/ccfang2/fdid?logo=GitHub)
![GitHub R package version](https://img.shields.io/github/r-package/v/ccfang2/fdid?logo=R)
![GitHub Repo stars](https://img.shields.io/github/stars/ccfang2/fdid?style=social)


> The R package `fdid` allows users to implement the method proposed in Fang and Liebl (2025)[^1]. In this paper, we present a novel functional perspective on Difference-in-Differences (DiD) that allows for honest inference using event study plots under violations of parallel trends and/or no-anticipation assumptions. We use the algorithm of so-called fast and fair simultaneous confidence band from Liebl and Reimherr (2023)[^2] to construct the simultaneous band in our plots.

## Installation

You can install the development version of `fdid` from [GitHub](https://github.com/) with:
      
``` r
# install.packages("devtools")
devtools::install_github("ccfang2/fdid")
```

## Example 1: Classic (Non-Honest) Simultaneous Inference

We hereby use event study estimates from Gallagher (2014)[^3] as an example. 

``` r
library(fdid)
data(Gdata)
fdid_scb_est <- fdid_scb(beta=Gdata$beta, cov=Gdata$cov, t0=Gdata$t0)
plot(fdid_scb_est, pos.legend="bottom", scale.legend=1.4)
```

We first load the data from our package, and then we use the function `fdid_scb()` to compute simultaneous confidence band from using the estimates of event study coefficients, covariances and reference time point, and finally we use the generic function `plot()` to derive a figure as follows.

![Example 1](man/figures/plot_scb.png)

As is seen, the reference time in this example is at event time -1. The gray area indicates the time span over which the treatment effect is uniformly significant under our simultaneous confidence band. Traditional confidence intervals are also plotted. The intervals fail to take into account the multiple testing problem, so they are typically narrower than our band. 

> **Note:** If you do not have estimates of event study coefficients and covariances, you may use our function `fdid()` to estimate them from using your original data. Our function allows the estimation under both non-staggered and staggered DiD designs. In particularly, we consider the negative weighting problem of estimating event study coefficients under staggered designs and use carefully chosen non-negative weights to sum up estimates from different treatment subgroups.

## Example 2: Honest Inference under Treatment anticipation

Following Example 1, we now suppose that, right after event time -3, there is an anticipation of treatment.

``` r
plot(fdid_scb_est, ta.t0=-3, pos.legend="bottom", scale.legend=1.4)
```

![Example 2](man/figures/plot_scb_ta.png)

With an anticipation right after event time -3, one may see the time span over which the treatment effect is uniformly significant shrinks. Under treatment anticipation, the null hypothesis for our test includes all values within the pink bounds.

## Example 3: Honest Inference under Differential Trend of Functional Relative Magnitudes

Following Example 1, we now suppose that, there is differential trend of functional relative magnitudes with control parameter $\overline{M}=1$.

``` r
plot(fdid_scb_est, frm.mbar=1, pos.legend="bottom", scale.legend=1.4)
```

![Example 3](man/figures/plot_scb_frm.png)

Under such a differential trend, the time span over which the treatment effect is uniformly significant is even smaller. It can be seen that the average absolute derivative of functional DiD estimate over pre-treatment time is quite large, so our bounds for honest inference is also large.

## Example 4: Honest Inference under Differential Trend of Functional Trend Restrictions

Following Example 1, we now suppose that, there is differential trend of functional trend restrictions with control parameter $M=2$.

``` r
plot(fdid_scb_est, ftr.m=2, pos.legend="bottom", scale.legend=1.4)
```

![Example 4](man/figures/plot_scb_ftr.png)

With this specific differential trend, the time span over which the treatment effect is uniformly significant is not so different from that without. It can be seen that the pre-trend, measured by the average of derivative of functional DiD estimate over pre-treatment time, is quite small. Hence, the inference result does not change too much from Example 1, even though the control parameter is not infinitesimal.

## Contact
Chencheng Fang, Email: [ccfang[at]uni-bonn.de](mailto:ccfang@uni-bonn.de),
Institute of Finance and Statistics, University of Bonn

## Reference
[^1]: Fang, C. and Liebl, D. (2025). Honest Difference-in-Differences Using Event Study Plots: A Functional Data Approach. Working Paper
[^2]: Liebl, D. and M. Reimherr (2023). Fast and fair simultaneous confidence bands for functional parameters. Journal of the Royal Statistical Society Series B: Statistical Methodology 85(3), 842–868
[^3]: Gallagher, J. (2014). Learning about an Infrequent Event: Evidence from Flood Insurance Take-Up in the United States. American Economic Journal: Applied Economics 6(3), 206–33.
