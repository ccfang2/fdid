#' Compute simultaneous confidence band for event study coefficients in functional DiD framework
#' @description The \code{fdid_scb} function is used to compute simultaneous confidence band for event study coefficients in a functional
#' framework, together with the point-wise confidence intervals at observable event time. We adopt the algorithm of the so-called fast and
#' fair simultaneous confidence bands. See Liebl and Reimherr (2023) for details.
#'
#' @param object an object of S3 class \code{"fdid"}, which is an output of the function \link{fdid}. If it is provided, arguments 'beta',
#' 'cov' and 't0' are ignored.
#' @param beta a numeric matrix with two columns. The first column is the estimates of event study coefficients and the second column is the
#' corresponding event time. The estimate of coefficient at reference time should be normalized to 0.
#' @param cov a numeric matrix of covariance estimates.
#' @param t0 a numeric scalar that indicates the reference time.
#' @param df degree of freedom for the t-distributed based bands. If NULL, Gaussian distributed bands are computed.
#' @param ci.alpha significance level for the pointwise confidence intervals.
#' @param scb.pre.alpha significance level for simultaneous confidence band in pre-treatment periods.
#' @param scb.post.alpha significance level for simultaneous confidence band in post-treatment periods.
#' @return The \code{fdid_scb} function returns a list containing three objects:
#' \describe{
#'  \item{scb}{a list that includes the functional estimate of event study coefficients and the simultaneous confidence band.}
#'  \item{ci}{a matrix that includes the point-wise estimates of event study coefficients and their confidence intervals.}
#'  \item{data}{a list that stores the original inputs of beta and cov.}
#' }
#' The output is an object of S3 class \code{"fdid_scb"}.
#' @import ffscb
#' @export
#'
#' @references Fang, C. and Liebl, D. (2025). Making Event Study Plots Honest: A Functional Data Approach to Causal Inference.
#' @references Liebl, D. and M. Reimherr (2023). Fast and fair simultaneous confidence bands for functional parameters. Journal of
#' the Royal Statistical Society Series B: Statistical Methodology 85(3), 842–868
#'
#' @seealso \link{fdid}, \link[ffscb]{confidence_band}, \link[ffscb]{confidence_band_fragm}, \link{cov_spline}
#'
#' @examples
#' data(LWdata)
#' fdid_scb_est1 <- fdid_scb(beta=LWdata$beta, cov=LWdata$cov, t0=LWdata$t0)
#'
#' data(simulated_stagger_example)
#' fdid_est2 <- fdid(data=simulated_stagger_data, treatment=simulated_stagger_treatment)
#' fdid_scb_est2 <- fdid_scb(object=fdid_est2)
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

  # compute spline interpolation
  len_t      <- length(betahat)
  len_spline <- 5*len_t

  timeVec_spline <- seq(min(timeVec), max(timeVec), length.out=len_spline)
  while (t0 %in% timeVec_spline) {len_spline <- len_spline+1; timeVec_spline <- seq(min(timeVec), max(timeVec), length.out=len_spline) }

  betahat_spline    <- spline(x=timeVec, y=betahat, n=len_spline, method="natural")$y
  betahat_splinefun <- splinefun(x=timeVec, y=betahat, method="natural")
  covhat_spline     <- cov_spline(cov=covhat, grid=timeVec, n_intrpl=len_spline)

  # if (any(is.na(covhat))) {
  #   covhat_spline <- apply(covhat, 2, function(col) { obs_idx <- !is.na(col); obs_event_t <- timeVec[obs_idx]; obs_col <- col[obs_idx]; spline_fit <- spline(obs_event_t, obs_col, xout = timeVec_spline, method = "natural")$y; spline_fit[timeVec_spline > max(obs_event_t) | timeVec_spline < min(obs_event_t)] <- NA; return(spline_fit)})
  #   covhat_spline <- apply(covhat_spline, 1, function(row) { obs_idx <- !is.na(row); obs_event_t <- timeVec[obs_idx]; obs_row <- row[obs_idx]; spline_fit <- spline(obs_event_t, obs_row, xout = timeVec_spline, method = "natural")$y; spline_fit[timeVec_spline > max(obs_event_t) | timeVec_spline < min(obs_event_t)] <- NA; return(spline_fit)})
  # } else{
  #   covhat_spline <- apply(covhat, 2, function(col) spline(x=timeVec, y=col, n=len_spline, method="natural")$y)
  #   covhat_spline <- apply(covhat_spline, 1, function(row) spline(x=timeVec, y=row, n=len_spline, method="natural")$y)
  # }

  # compute tau function for pre-treatment periods
  betahat_spline_pre     <- betahat_spline[which(timeVec_spline<=t0)]
  covhat_spline_pre      <- covhat_spline[which(timeVec_spline<=t0), which(timeVec_spline<=t0)]
  diag_covhat_spline_pre <- diag(covhat_spline_pre)
  hat.tau_pre            <- ffscb::cov2tau_fun(covhat_spline_pre)

  # compute simultaneous confidence band for pre-treatment periods
  if(is.null(df)){
    if (any(is.na(covhat_spline_pre))) {
      scb_pre <- ffscb::confidence_band_fragm(x=betahat_spline_pre, diag.cov.x=diag_covhat_spline_pre, tau=hat.tau_pre, type = "FFSCB.z", conf.level=1-2*scb.pre.alpha, n_int=1)
    } else {
      scb_pre <- ffscb::confidence_band(x=betahat_spline_pre, cov.x=covhat_spline_pre, tau=hat.tau_pre, type="FFSCB.z", conf.level=1-2*scb.pre.alpha, n_int=1)
    }
  } else{
    if (any(is.na(covhat_spline_pre))) {
      scb_pre <- ffscb::confidence_band_fragm(x=betahat_spline_pre, diag.cov.x=diag_covhat_spline_pre, tau=hat.tau_pre, df=df, type = "FFSCB.t", conf.level=1-2*scb.pre.alpha, n_int=1)
    } else {
      scb_pre <- ffscb::confidence_band(x=betahat_spline_pre, cov.x=covhat_spline_pre, tau=hat.tau_pre,  df=df, type="FFSCB.t", conf.level=1-2*scb.pre.alpha, n_int=1)
    }
  }

  # compute tau function for post-treatment periods
  betahat_spline_post     <- betahat_spline[which(timeVec_spline>=t0)]
  covhat_spline_post      <- covhat_spline[which(timeVec_spline>=t0), which(timeVec_spline>=t0)]
  diag_covhat_spline_post <- diag(covhat_spline_post)
  hat.tau_post            <- ffscb::cov2tau_fun(covhat_spline_post)

  # compute simultaneous confidence band for post-treatment periods
  if(is.null(df)){
    if (any(is.na(covhat_spline_post))) {
      scb_post <- ffscb::confidence_band_fragm(x=betahat_spline_post, diag.cov.x=diag_covhat_spline_post, tau=hat.tau_post, type = "FFSCB.z", conf.level=1-scb.post.alpha, n_int=1)
    } else {
      scb_post <- ffscb::confidence_band(x=betahat_spline_post, cov.x=covhat_spline_post, tau=hat.tau_post, type="FFSCB.z", conf.level=1-scb.post.alpha, n_int=1)
    }
  } else{
    if (any(is.na(covhat_spline_post))) {
      scb_post <- ffscb::confidence_band_fragm(x=betahat_spline_post, diag.cov.x=diag_covhat_spline_post, tau=hat.tau_post, df=df, type = "FFSCB.t", conf.level=1-scb.post.alpha, n_int=1)
    } else {
      scb_post <- ffscb::confidence_band(x=betahat_spline_post, cov.x=covhat_spline_post, tau=hat.tau_post,  df=df, type="FFSCB.t", conf.level=1-scb.post.alpha, n_int=1)
    }
  }

  idx_pos          <- which(timeVec_spline>=t0)[1]
  scb              <- rbind(scb_pre, 0, scb_post)
  x                <- c(timeVec_spline[1:(idx_pos-1)],t0,timeVec_spline[idx_pos:len_spline])

  scb_ub_splinefun <- splinefun(x=x, y=scb[,(ncol(scb)-1)], method="natural")
  scb_lb_splinefun <- splinefun(x=x, y=scb[,ncol(scb)], method="natural")

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

