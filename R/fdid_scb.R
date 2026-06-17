#' Compute simultaneous confidence bands for event study coefficients in a functional DiD framework
#' @description The \code{fdid_scb} function is used to compute simultaneous confidence bands for event study coefficients in a functional
#' framework. Specifically, we compute an infimum-based simultaneous confidence band in pre-treatment period by parametric bootstrap, and a supremum-based
#' simultaneous confidence band in post-treatment period by an algorithm of Kac-Rice formula proposed in \href{https://academic.oup.com/jrsssb/article/85/3/842/7133768}{Liebl and Reimherr (2023)}.
#'
#' @param object an object of S3 class \code{"fdid"}, which is an output of the function \link{fdid}. If it is provided, arguments 'beta',
#' 'cov' and 't0' are ignored.
#' @param beta a numeric matrix with two columns. The first column is the estimates of event study coefficients and the second column is the
#' corresponding event time. The estimate of coefficient at reference time should be normalized to 0.
#' @param cov a numeric matrix of covariance estimates.
#' @param t0 a numeric scalar that indicates the reference time.
#' @param df degree of freedom for the t-distributed based band in post-treatment periods. If NULL, Gaussian distributed bands are computed.
#' @param ci.alpha significance level for the point-wise confidence intervals.
#' @param scb.pre.alpha significance level for the infimum-based simultaneous confidence band in pre-treatment periods.
#' @param scb.post.alpha significance level for the supremum-based simultaneous confidence band in post-treatment periods.
#' @return The \code{fdid_scb} function returns a list containing three objects:
#' \describe{
#'  \item{scb}{a list that includes the functional estimate of event study coefficients and the simultaneous confidence bands.}
#'  \item{ci}{a matrix that includes the point-wise estimates of event study coefficients and their confidence intervals.}
#'  \item{data}{a list that stores the original inputs of beta and cov.}
#' }
#' The output is an object of S3 class \code{"fdid_scb"}.
#' @import ffscb
#' @export
#'
#' @references Fang, C. and Liebl, D. (2026). Making Event Study Plots Honest: A Functional Data Approach to Causal Inference. \href{https://arxiv.org/abs/2512.06804}{arXiv:2512.06804}.
#' @references Liebl, D. and M. Reimherr (2023). Fast and fair simultaneous confidence bands for functional parameters. Journal of
#' the Royal Statistical Society Series B: Statistical Methodology 85(3), 842–868
#'
#' @seealso \link{fdid}, \link[ffscb]{confidence_band}, \link[ffscb]{confidence_band_fragm}, \link{cov_spline}
#'
#' @examples
#' data(LWdata)
#' fdid_scb_est <- fdid_scb(beta=LWdata$beta, cov=LWdata$cov, t0=LWdata$t0)
fdid_scb <- function(object=NULL,
                     beta=NULL,
                     cov=NULL,
                     t0=NULL,
                     df=NULL,
                     ci.alpha=0.05,
                     scb.pre.alpha=0.05,
                     scb.post.alpha=0.05) {

  # check conditions
  if (!is.null(object) && !base::inherits(object,"fdid")) stop("The input 'object' should be either NULL or an output of function 'fdid'.")
  if (!is.null(object)) warning("The input 'object' is provided, so 'beta', 'cov' and 't0' are ignored.")

  if (!is.null(beta) && !is.matrix(beta)) stop("The input 'beta' should be either NULL or a numeric matrix.")
  if (is.null(beta) && is.null(object)) stop("If any of 'beta', 'cov' and 't0' is NULL, 'object' must be provided.")
  if (is.matrix(beta) && !is.numeric(beta)) stop("If 'beta' is not NULL, it should be a numeric matrix.")
  if (is.matrix(beta) && any(is.na(beta))) stop("The input 'beta' should contain no missing values.")
  if (!is.null(cov) && !is.matrix(cov)) stop("The input 'cov' should be eitherNULL or a matrix.")
  if (is.null(cov) && is.null(object)) stop("If any of 'beta', 'cov' and 't0' is NULL, 'object' must be provided.")
  if (!is.null(t0) && (!is.numeric(t0) || length(t0) != 1)) stop("The input 't0' should be either NULL or a numeric scalar.")
  if (!is.null(df) && (!is.numeric(df) || length(df) != 1 || df <= 0)) stop("The input 'df' should be either NULL or a positive integer scalar.")

  # extract data
  if(is.null(object)) {

    # check further conditions
    if (!(t0 %in% beta[,2])) stop("t0 should be in the time vector.")
    if (beta[,1][which(beta[,2]==t0)]!=0) stop("beta at t0 should be normalized to 0.")
    if (any(cov[which(beta[,2]==t0),]!=0) && any(covhat[,which(timeVec==t0)]!=0)) stop("Rows and columns of cov at t0 should all be normalized to 0.")

    #beta[,2]       <- beta[,2]-t0
    t_order        <- order(beta[,2])
    beta           <- beta[t_order,]
    colnames(beta) <- c(colnames(beta)[1], "event_t")

    cov            <- cov[t_order, t_order]
    colnames(cov)  <- rownames(cov) <- beta[,"event_t"]

    betahat        <- beta[,1]
    covhat         <- cov
    timeVec        <- beta[,2]
    t0             <- t0
    df             <- df

  } else {

    beta           <- object$beta$coef
    cov            <- object$beta$cov

    betahat        <- beta[,1]
    covhat         <- cov
    timeVec        <- beta[,2]
    t0             <- object$t0
    df             <- object$df
  }

  # estimate point-wise confidence interval
  if (is.null(df)) {
    ci_upper <- betahat+qnorm(p=1-ci.alpha/2)*sqrt(diag(covhat))
    ci_lower <- betahat-qnorm(p=1-ci.alpha/2)*sqrt(diag(covhat))
  } else {
    ci_upper <- betahat+qt(p=1-ci.alpha/2,df=df)*sqrt(diag(covhat))
    ci_lower <- betahat-qt(p=1-ci.alpha/2,df=df)*sqrt(diag(covhat))
  }

  # compute infimum-based simultaneous confidence band for pre-treatment periods via parametric bootstrapping
  timeVec_pre     <- timeVec[which(timeVec<=t0)]
  betahat_pre     <- betahat[which(timeVec<=t0)]
  covhat_pre      <- covhat[which(timeVec<=t0), which(timeVec<=t0)]
  diag_covhat_pre <- diag(covhat_pre)

  betahat_spline_pre     <- spline(x=timeVec_pre, y=betahat_pre, n=7*length(timeVec_pre), method="natural")$y
  betahat_splinefun_pre  <- splinefun(x=timeVec_pre, y=betahat_pre, method="natural")
  covhat_spline_pre      <- cov_spline(cov=covhat_pre, grid=timeVec_pre, n_intrpl=7*length(timeVec_pre))
  diag_covhat_spline_pre <- diag(covhat_spline_pre)

  B <- 300 #500
  T_boot <- numeric(B)

  set.seed(450)

  # for (b in 1:B) {
  #   betahat_pre_b <- MASS::mvrnorm(1, mu=betahat_pre, Sigma=covhat_pre)
  #   betahat_spline_pre_b <- spline(betahat_pre_b, n=len_spline_pre, method="natural")$y
  #   t_b <- (betahat_spline_pre_b - betahat_spline_pre)[-len_spline_pre] / sqrt(diag_covhat_spline_pre)[-len_spline_pre]
  #   T_boot[b] <- min(abs(t_b))
  # }

  for (b in 1:B) {
    betahat_pre_b <- MASS::mvrnorm(1, mu=betahat_pre, Sigma=covhat_pre)
    t_b <- (betahat_pre_b - betahat_pre)[-length(betahat_pre)] / sqrt(diag(covhat_pre))[-length(betahat_pre)]
    T_boot[b] <- min(abs(t_b))
  }

  c_alpha <- quantile(T_boot, probs=1-2*scb.pre.alpha, type=6)
  scb_pre <- cbind(betahat_spline_pre, betahat_spline_pre+c_alpha*sqrt(diag_covhat_spline_pre), betahat_spline_pre-c_alpha*sqrt(diag_covhat_spline_pre))

  # compute supremum-based simultaneous confidence band for post-treatment periods by Kac-Rice formula
  # leave the first interval after the reference period empty
  timeVec_post     <- timeVec[which(timeVec>=timeVec[which(timeVec==t0)+1])]
  betahat_post     <- betahat[which(timeVec>=timeVec[which(timeVec==t0)+1])]
  covhat_post      <- covhat[which(timeVec>=timeVec[which(timeVec==t0)+1]), which(timeVec>=timeVec[which(timeVec==t0)+1])]
  diag_covhat_post <- diag(covhat_post)

  betahat_spline_post     <- spline(x=timeVec_post, y=betahat_post, n=7*length(timeVec_post), method="natural")$y
  betahat_splinefun_post  <- splinefun(x=timeVec_post, y=betahat_post, method="natural")
  covhat_spline_post      <- cov_spline(cov=covhat_post, grid=timeVec_post, n_intrpl=7*length(timeVec_post))
  diag_covhat_spline_post <- diag(covhat_spline_post)
  hat.tau_post            <- ffscb::cov2tau_fun(covhat_spline_post) # there should be no zero in hat.tau_post

  if(is.null(df)){
      scb_post <- ffscb::confidence_band(x=betahat_spline_post, cov.x=covhat_spline_post, tau=hat.tau_post, type="FFSCB.z", conf.level=1-scb.post.alpha, n_int=1)
  } else{
      scb_post <- ffscb::confidence_band(x=betahat_spline_post, cov.x=covhat_spline_post, tau=hat.tau_post,  df=df, type="FFSCB.t", conf.level=1-scb.post.alpha, n_int=1)
  }

  # derive the spline functions
  scb_ub_splinefun_pre  <- splinefun(x=seq(timeVec_pre[1], t0, len=7*length(timeVec_pre)), y=scb_pre[,(ncol(scb_pre)-1)], method="natural")
  scb_lb_splinefun_pre  <- splinefun(x=seq(timeVec_pre[1], t0, len=7*length(timeVec_pre)), y=scb_pre[,ncol(scb_pre)], method="natural")

  scb_ub_splinefun_post  <- splinefun(x=seq(timeVec_post[1], timeVec_post[length(timeVec_post)], len=7*length(timeVec_post)), y=scb_post[,(ncol(scb_post)-1)], method="natural")
  scb_lb_splinefun_post  <- splinefun(x=seq(timeVec_post[1], timeVec_post[length(timeVec_post)], len=7*length(timeVec_post)), y=scb_post[,ncol(scb_post)], method="natural")

  scb_ub_splinefun <- function(t) ifelse(t<=t0, scb_ub_splinefun_pre(t), ifelse(t>= timeVec_post[1], scb_ub_splinefun_post(t), NA))
  scb_lb_splinefun <- function(t) ifelse(t<=t0, scb_lb_splinefun_pre(t), ifelse(t>= timeVec_post[1], scb_lb_splinefun_post(t), NA))
  betahat_splinefun <- function(t) (scb_ub_splinefun(t)+scb_lb_splinefun(t))/2

  # final output
  final_output <- list(scb=list(betahat= betahat_splinefun,
                                scb_ub= scb_ub_splinefun,
                                scb_lb= scb_lb_splinefun,
                                event_t=timeVec),
                       ci=cbind(betahat, ci_upper, ci_lower, event_t=timeVec),
                       data=list(beta=beta, cov=cov, t0=t0, df=df, ci.alpha=ci.alpha, scb.pre.alpha=scb.pre.alpha, scb.post.alpha=scb.post.alpha))

  class(final_output) <- "fdid_scb"
  return(final_output)

  options(warn=1)
}

# Extension: In staggered adoption design, we focus on a common time window.

