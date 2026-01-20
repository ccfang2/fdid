# fdid: Making Event Study Plots Honest 
## A Functional Data Approach to Causal Inference <img src="man/figures/fdid_badge.png" align="right" alt="" width="155" />

[![PRs Welcome](https://img.shields.io/badge/PRs-welcome-brightgreen.svg?style=flat-square)](https://makeapullrequest.com)
![GitHub last commit](https://img.shields.io/github/last-commit/ccfang2/fdid?logo=GitHub)
![GitHub repo size](https://img.shields.io/github/repo-size/ccfang2/fdid?logo=GitHub)
![GitHub R package version](https://img.shields.io/github/r-package/v/ccfang2/fdid?logo=R)
![GitHub Repo stars](https://img.shields.io/github/stars/ccfang2/fdid?style=social)


The R package `fdid` allows users to implement the method proposed in Fang and Liebl (2026)[^1]. In this paper, we present a novel functional perspective on Difference-in-Differences (DiD) that allows for honest inference using event study plots under violations of parallel trends and/or no-anticipation assumptions. Specifically, we compute an infimum-based simultaneous confidence band in the pre-treatment period by parametric bootstrap, and a supremum-based simultaneous confidence band in the post-treatment period by the algorithm of Kac-Rice formula proposed in Liebl and Reimherr (2023)[^2]. Additionally, by contrast to classical reference line in traditional event study plots, we derive an honest reference band, accounting for potential biases from the violation of parallel trends or no-anticipation assumption, when making inference.

By doing so, we turn traditional event study plots into rigorous honest causal inference tools through equivalence and relevance testing: Honest reference band can be validated via equivalence testing in the pre-anticipation period, and honest causal effects can be tested through relevance testing using the honest reference band in the post-treatment period.

> You may find a presentation of an early version of the paper in my YouTube [video](https://www.youtube.com/watch?v=h0KCv8y9Apw). Also, users can adjust their honest reference bands interactively via our [Shiny app](https://ccfang2.shinyapps.io/fdidHonestInference/).

## Installation

You can install the development version of `fdid` from [GitHub](https://github.com/) with:
      
``` r
# install.packages("devtools")
devtools::install_github("ccfang2/fdid")
```

## Classical Event Study Plot

We hereby use event study estimates from Gallagher (2014)[^3]. The following is the traditional event study plot displaying pointwise 95% confidence intervals.

```r
library(fdid)
data(Gdata)
Gdata$beta[,"event_t"] <- Gdata$beta[,"event_t"]- Gdata$t0 #Recenter the event time on 0

fdid_scb_est <- fdid_scb(beta=Gdata$beta, cov=Gdata$cov, t0=0)
par(cex.axis = 1.4, cex.lab = 1.4, cex.main = 1.4, family="Times")
EventStudyPlot_Classical(fdid_scb_est, pos.legend="bottom", scale.legend=1.4)
```

<p align="center">
<img src="man/figures/plot_ci.png" width="80%">
</p>

> In the original dataset `Gdata` from Gallagher (2014)[^3], they already consider a potential anticipation starting after event time -1, which is used as the reference time point. However, in our approach, we suggest always using event time 0 as the reference time point, and then derive the honest reference band that considers potential violation of no-anticipation or parallel trends assumption. To accommodate our approach, we thus recenter the event time on 0 in the dataset `Gdata`.

> The function `fdid_scb()` is used to compute simultaneous confidence bands from using the estimates of event study coefficients, covariances and reference time point, which will be used in the honest causal inference. 

However, classical event study plots---such as the one above---suffer from at least three important limitations. 

 1. First, they typically display pointwise confidence intervals that do not account for multiple testing across event times.
 2. Second, and of particular practical importance, they give the impression that the parallel trends and no-anticipation assumptions can be validated in the pre-treatment period when showing insignificant pre-treatment estimates. However, this is a classical argument from ignorance as the failure to reject the null hypothesis of no pre-treatment effects (i.e. absence of differences in time trends or anticipatory effects) does not imply that the parallel trends and no-anticipation as- sumptions hold.
 3. Third, and equally important from a practical standpoint, honest inference methods—such as those developed by Rambachan and Roth (2023)[^4]—cannot be integrated into standard event study plots, limiting their usefulness for credible causal inference.

We address these limitations by introducing a functional-data perspective on DiD. The key idea is to model the underlying time-series processes in continuous time—an assumption already implicit in many empirical DiD studies, where pointwise event-study estimates and confidence intervals are connected by straight lines across event times. Our estimator builds directly on standard panel-data structures and is pointwise identical to the classical panel estimator, making the approach straightforward to implement in empirical applications. We allow for both additional control variables and staggered treatment adoption.

1. First, the Gaussian process result provides the foundation for constructing simultaneous confidence bands for the DiD parameter (i.e., the event-study coefficients) across the full continuum of event times. Compared to conventional pointwise inference, our approach offers a powerful and more credible alternative, explicitly accounting for the multiple-testing problem inherent in event-study analyses.
2. Second, our infimum-based simultaneous confidence bands enable formal validation of honest reference bands in the pre-treatment period via equivalence testing, allowing researchers to assess the plausibility of the parallel trends and no-anticipation assumptions in a rigorous statistical manner. 
3. Third, our supremum-based simultaneous confidence bands support honest inference in the post-treatment period through relevance testing, directly integrating existing approaches to honest DiD inference into the event-study framework. 


## Simultaneous Confidence Bands

We can therefore transform the traditional event study plot into a rigorous honest inference tool with the infimum-based 90% simultaneous confidence band in pre-treatment period and supremum-based 95% simultaneous confidence band in post-treatment period. The infimum-based 90% simultaneous confidence band is for performing equivalence testing at significance level 5% (see Section 3.3 in Fang and Liebl (2026)[^1]), i.e. validating the honest reference band in the pre-anticipation period; and the supremum-based 95% simultaneous confidence band is for performing relevance testing at significance level 5% (see Section 3.1 in Fang and Liebl (2026)[^1]), i.e. uniformly and honestly testing causal inference in the post-treatment period. The following is the new plot using simultaneous confidence bands.

``` r
par(cex.axis = 1.4, cex.lab = 1.4, cex.main = 1.4, family="Times")
plot(fdid_scb_est, pos.legend="bottom", scale.legend=1.4, note.pre=FALSE, ci.post=TRUE)
```

<p align="center">
<img src="man/figures/plot_scb.png" width="80%">
</p>

We use the generic function `plot()` to derive the plot. In the post-treatment period, the supremum-based 95% simultaneous confidence band is wider than the classical 95% confidence intervals, because the pointwise intervals fail to take into account the multiple testing. The treatment effect is uniformly significant in the simultaneous causal inference using the classical reference line over event time [0, 9.7]. 

> We do not perform validation in the pre-treatment period, because we only use the classical reference line in the simultaneous inference above.

> If you do not have estimates of event study coefficients and covariances, you may use our function `fdid()` to estimate them from using your original data. Our function allows the estimation under both non-staggered and staggered DiD designs. In particularly, we consider the negative weighting problem of estimating event study coefficients under staggered designs and use carefully chosen non-negative weights to sum up estimates from different treatment subgroups.

To conduct honest inference using the plot above, we need to derive the honest reference band under violations of identification assumptions. Below are two examples of deriving honest reference band under the violation of no-anticipation and parallel trends assumption respectively.

## Example 1: Honest Reference Band under Violation of No-anticipation Assumption

We now suppose that, after event time -2, there is an anticipation of treatment. We use control parameters $S_{u}=1.4$ and $S_{\ell}=2.3$ to derive the reference band (see equation (36) in Fang and Liebl (2026)[^1] for details on the control parameters).

``` r
par(cex.axis = 1.4, cex.lab = 1.4, cex.main = 1.4, family="Times")
plot(fdid_scb_est, ta.ts=-2, ta.s=c(1.4,2.3), pos.legend="bottom", scale.legend=1.4, ci.post=TRUE, ref.band.pre = TRUE)
```

<p align="center">
<img src="man/figures/plot_scb_ta.png" width="80%">
</p>

With an anticipation after event time -2, one may see that the treatment effect is still uniformly significant over event time [0.4, 9,0]. With the given control parameters, the reference band can be validated at the significance level 5%, since the infimum-based 90% simultaneous confidence band strictly lies within the reference band in the pre-anticipation period (see Section 3.3 in Fang and Liebl (2026)[^1] for details). The result shows that the treatment effect in Gallagher(2014)[^3] is robust under the considered treatment anticipation.

## Example 2: Honest Reference Band under Violation of Parallel Trends Assumption

We now suppose that, there is differential trend. We use control parameters $M_{u}=0.25$ and $M_{\ell}=0.25$ to derive the reference band (see equation (37) in Fang and Liebl (2026)[^1] for details on the control parameters).

``` r
par(cex.axis = 1.4, cex.lab = 1.4, cex.main = 1.4, family="Times")
plot(fdid_scb_est, frmtr.m=c(0.25,0.25), pos.legend="bottom", scale.legend=1.4, ci.post=TRUE, ref.band.pre = TRUE)
```

<p align="center">
<img src="man/figures/plot_scb_frmtr.png" width="80%">
</p>

In the plot above, although the reference band cannot be validated at the significance level 5% due to high data variability, it captures the visible upward pre-trend with a width comparable to that of the infimum-based band, providing substantive justification. Using this reference band, we find that the treatment effect is still uniformly significant over event time [0, 7.7]. The result shows that the treatment effect in Gallagher(2014)[^3] is robust under the considered violation of parallel trends assumption.

> In some cases, validating a given reference band can be challenging, as doing so may require selecting a very wide reference band—thereby making subsequent testing in the post-treatment period overly conservative. Such non-rejection of the equivalence null hypothesis (see Section 3.3 in Fang and Liebl (2026)[^1] for details) often reflects limited sample size or high variability, and must be viewed as a lack of evidence against the null, not confirmation of it. Thus, a reference band failing to pass the validation can still be used for honest inference when its specification can be supported by domain-specific justification.

## Contact
Chencheng Fang, Email: [ccfang[at]uni-bonn.de](mailto:ccfang@uni-bonn.de), Hausdorff Center for Mathematics; Institute of Finance and Statistics, University of Bonn


[^1]: Fang, C. and Liebl, D. (2026). Making Event Study Plots Honest: A Functional Data Approach to Causal Inference. [arXiv:2512.06804](https://arxiv.org/abs/2512.06804).
[^2]: Liebl, D. and M. Reimherr (2023). Fast and fair simultaneous confidence bands for functional parameters. Journal of the Royal Statistical Society Series B: Statistical Methodology 85(3), 842–868.
[^3]: Gallagher, J. (2014). Learning about an Infrequent Event: Evidence from Flood Insurance Take-Up in the United States. American Economic Journal: Applied Economics 6(3), 206–33.
[^4]: Rambachan, A. and J. Roth (2023). A more credible approach to parallel trends. The Review
of Economic Studies 90 (5), 2555–2591.
