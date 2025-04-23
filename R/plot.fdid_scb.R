#' Custom Plot Method for objects of class \code{"fdid_scb"}
#' @description This function is an S3 method for \link[base]{plot}, specifically designed for objects of class \code{"fdid_scb"}.
#' It plots simultaneous confidence band with bounds for honest inference.
#'
#' @param object an object of class \code{"fdid_scb"}. The object to be plotted.
#' @param ta.t0 a numeric value that indicates the time, right after which there is an anticipation of treatment. It should be
#' greater than the minimal event time and no greater than t0 in \code{object}. If it is NULL, there is no anticipation.
#' @param ta.s a numeric value of control parameter for honest inference under treatment anticipation. If NULL, the critical value \eqn{t_{1-\alpha/2, df}} or \eqn{z_{1-\alpha/2}} is applied with the same \eqn{\alpha} used for constructing simultaneous confidence band. If no treatment anticipation is defined in \code{ta.t0}, the value of \code{ta.s} is ignored.
#' @param frm.mbar a numeric value of control parameter for honest inference under functional relative magnitudes.
#' @param ftr.m a numeric value of control parameter for honest inference under functional trend restrictions.
#' @param frmtr.mbar a numeric value of control parameter for combining functional relative magnitudes and functional trend restrictions.
#' @param post.trt a logical value. If TRUE, only bounds for post-treatment periods are plotted and used for honest inference.
#' @param ci a logical value. If TRUE, the point-wise confidence intervals are also plotted.
#' @param pos.legend a character value of "top" or "bottom" that indicates the position of legend. If NULL, the legend is not printed.
#' @param ... Additional arguments to be passed to \link[base]{plot}.
#'
#' @return The function returns a plot with simultaneous confidence band for event study coefficients in a functional framework,
#' together with bounds for honest inference, if properly defined.
#' @note The inference result around the reference time point (i.e. between two event time closest to the reference time) should be treated with caution.
#'
#' @import pracma
#' @export
#'
#' @references Fang, C. and Liebl, D. (2025). Honest Causal Inference with Difference-in-Differences: A Functional Data Perspective.
#' @seealso \link{fdid_scb}
#'
#' @examples
#' data(LWdata)
#' fdid_scb_est <- fdid_scb(beta=LWdata$beta, cov=LWdata$cov, t0=LWdata$t0)
#' cat("The reference time is ", LWdata$t0, ". If not NULL, the input 'ta.t0' in function 'plot' should be smaller than this value.", sep="")
#'
#' ## simultaneous inference
#' plot(fdid_scb_est)
#'
#' ## honest inference under treatment anticipation
#' plot(fdid_scb_est, ta.t0=-4)
#'
#' ## honest inference under violation of parallel trends
#' plot(fdid_scb_est, frm.mbar=0.3)
plot.fdid_scb <- function(object,
                          ta.t0=NULL,
                          ta.s=NULL,
                          frm.mbar=NULL,
                          ftr.m=NULL,
                          frmtr.mbar=NULL,
                          post.trt=FALSE,
                          ci=TRUE,
                          pos.legend="top",
                          ...) {

  # give out a warn
  warning_note <- paste("The inference result around the reference time t=", object$data$t0, "(i.e. between two event time closest to the reference time) should be treated with caution.")
  warning(warning_note)

  options(warn=-1)

  # check conditions
  if (!base::inherits(object,"fdid_scb")) stop("The input 'object' should be an output of function 'fdid_scb'.")
  if (!is.null(ta.t0) && (!is.numeric(ta.t0) || length(ta.t0) != 1)) stop("The input 'ta.t0' should be either NULL or a numeric scalar.")
  if (!is.null(ta.t0) && !(ta.t0 %in% object$scb$event_t[which(object$scb$event_t<=object$data$t0)])) stop("If not NULL, the input 'ta.t0' should be among the pre-treatment event time in 'object'.")
  if (!is.null(ta.t0) && ta.t0 == object$scb$event_t[1] && !all(sapply(list(frm.mbar, ftr.m, frmtr.mbar), is.null))) stop("If 'ta.t0' is defined to be the first event time, 'frm.mbar', 'ftr.m' and 'frmtr.mbar' must be NULL, because there is no available data for computing pre-trend differences.")
  if (!is.null(ta.s) && (!is.numeric(ta.s) || length(ta.s) != 1 || ta.s < 0)) stop("The input 'ta.s' should be either NULL or a numeric non-negative scalar.")
  if (!is.null(frm.mbar) && (!is.numeric(frm.mbar) || length(frm.mbar) != 1 || frm.mbar < 0)) stop("The input 'frm.mbar' should be either NULL or a numeric non-negative scalar.")
  if (!is.null(ftr.m) && (!is.numeric(ftr.m) || length(ftr.m) != 1 || ftr.m < 0)) stop("The input 'ftr.m' should be either NULL or a numeric non-negative scalar.")
  if (!is.null(frmtr.mbar) && (!is.numeric(frmtr.mbar) || length(frmtr.mbar) != 1 || frmtr.mbar < 0)) stop("The input 'frmtr.mbar' should be either NULL or a numeric non-negative scalar.")
  if (!is.null(frm.mbar) && !is.null(ftr.m) && !is.null(frmtr.mbar)) stop("Exactly one of inputs 'frm.mbar', 'ftr.m' and 'frmtr.mbar' can have a value, or all three are NULL.")
  if (!is.logical(post.trt)) stop("The input 'post.trt' should be logical.")
  if (!is.logical(ci)) stop("The input 'ci' should be logical.")
  if (!is.null(pos.legend) && !pos.legend %in% c("top", "bottom")) stop("The input 'pos.legend' must be 'top', 'bottom', or NULL.")

  # extract data from object
  betahat <- object$data$beta[,1]
  covhat  <- object$data$cov
  timeVec <- object$data$beta[,2]

  betahat_splinefun <- object$scb$betahat
  scb_ub_splinefun  <- object$scb$scb_ub
  scb_lb_splinefun  <- object$scb$scb_lb

  ci_upper <- object$ci[,"ci_upper"]
  ci_lower <- object$ci[,"ci_lower"]

  start                     <- min(object$data$beta[,2])
  end                       <- max(object$data$beta[,2])
  t0                        <- object$data$t0
  if (is.null(ta.t0)) ta.t0 <- t0

  # regular statistical inference
  if (ta.t0==t0 & is.null(frm.mbar) & is.null(ftr.m) & is.null(frmtr.mbar)) {

    # find the time spans of statistical significance
    n.int <- 200
    roots <- vector("numeric", length = n.int)
    for (i in 1:n.int) {
      fun      <- function(x) scb_ub_splinefun(x)*scb_lb_splinefun(x)
      interval <- c(start+(i-1)*(end-start)/n.int, start+(i)*(end-start)/n.int)
      root     <- try(uniroot(fun, interval = interval)$root, silent = TRUE)
      roots[i] <- if(inherits(root, "try-error")) {NA} else root
    }
    #roots <- roots[!is.na(roots)]
    roots <- sort(unique(roots[!is.na(roots)]))

    roots_vec    <- sort(unique(c(start,roots,end)))
    stat_sig_vec <- rep(NA, length(roots_vec)-1)
    for (i in 1:(length(roots_vec)-1)) {
      test_point      <- (roots_vec[i]+roots_vec[i+1])/2
      scb_ub_sign     <- sign(scb_ub_splinefun( test_point))
      scb_lb_sign     <- sign(scb_lb_splinefun( test_point))
      stat_sig_vec[i] <- (scb_ub_sign+scb_lb_sign)/2
    }

    x_vals  <- seq(start, end, length.out = 5*length(timeVec))
    y_range <- range(c(betahat_splinefun(x_vals), scb_ub_splinefun(x_vals), scb_lb_splinefun(x_vals)), na.rm=TRUE)

    # draw plot
    plot(timeVec, betahat,  type = "p", ylim =y_range, xlab = "Event Time", ylab = "", pch = 16, col = "blue", cex.axis=1.4, cex.lab=1.4, ...)
    title(ylab=expression(hat(beta)), mgp=c(2.1,1,0), cex.lab=1.4)
    if(isTRUE(ci)) arrows(timeVec, ci_lower, timeVec, ci_upper, angle = 90, code = 3, length = 0.025, col = "blue4")
    grid()

    curve(betahat_splinefun, add=TRUE, lty=1)
    curve(scb_ub_splinefun, add=TRUE, lty=1)
    curve(scb_lb_splinefun, add=TRUE, lty=1)

    if(isTRUE(post.trt)) {
      segments(x0=t0, y0=0, x1=end, y1=0, lty=3, lwd=4, col="red")

      if(!(t0 %in% roots_vec)) {
        roots_vec <- sort(c(roots_vec, t0))
        stat_sig_vec <- append(stat_sig_vec, stat_sig_vec[which(roots_vec==t0)-1], after=which(roots_vec==t0)-1-1)
      }

      for (i in 1:(length(roots_vec)-1)) {
        if (stat_sig_vec[i]==1 & roots_vec[i]>=t0 ) {
          rect(roots_vec[i], y_range[1], roots_vec[i+1], y_range[2], col=rgb(0.5, 0.5, 0.5, alpha=0.7), border=NA)
        } else {
          if (stat_sig_vec[i]==-1 & roots_vec[i]>=t0) {
            rect(roots_vec[i], y_range[1], roots_vec[i+1], y_range[2], col=rgb(0.5, 0.5, 0.5, alpha=0.4), border=NA)
          }
        }
      }
    } else {
      segments(x0=start, y0=0, x1=end, y1=0, lty=3, lwd=4, col="red")
      for (i in 1:(length(roots_vec)-1)) {
        if (stat_sig_vec[i]==1) {
          rect(roots_vec[i], y_range[1], roots_vec[i+1], y_range[2], col=rgb(0.5, 0.5, 0.5, alpha=0.7), border=NA)
        } else {
          if (stat_sig_vec[i]==-1) {
            rect(roots_vec[i], y_range[1], roots_vec[i+1], y_range[2], col=rgb(0.5, 0.5, 0.5, alpha=0.4), border=NA)
          }
        }
      }
    }

    # write a legend for the plot
    if(!is.null(pos.legend)) {

      xlim <- par("usr")[1:2]
      ylim <- par("usr")[3:4]

      text              <- "Reference value          "
      text_width        <- strwidth(text)
      fixed_line_width  <- strwidth("Reference")
      extra_padding     <- strwidth("   ", cex=1)
      rect_width        <- text_width+extra_padding+fixed_line_width
      rect_height       <- diff(ylim)*0.08
      rect_x1           <- mean(xlim)-rect_width/2
      rect_x2           <- mean(xlim)+rect_width/2

      if(pos.legend=="top")    {rect_y2 <- ylim[2]-rect_height*0.5; rect_y1 <- rect_y2-rect_height}
      if(pos.legend=="bottom") {rect_y2 <- ylim[1]+rect_height*0.5; rect_y1 <- rect_y2+rect_height}

      rect(rect_x1, rect_y1, rect_x2, rect_y2, col = "white", border = "black")

      left_padding <- rect_width*0.05
      line_x1      <- rect_x1+left_padding
      line_x2      <- rect_x1+left_padding+fixed_line_width

      if(pos.legend=="top")    {line_y <- rect_y1+rect_height*0.4}
      if(pos.legend=="bottom") {line_y <- rect_y1-rect_height*0.4}

      segments(line_x1, line_y, line_x2, line_y, col = "red", lwd = 4, lty = 3)
      text(line_x2 + left_padding, line_y, text, cex = 1, adj = 0)
    }
  } else {

    # honest inference
    if(is.null(frm.mbar) & is.null(ftr.m) & is.null(frmtr.mbar)) {

      if(is.null(ta.s)) {
        ta_ub <- ci_upper[names(ci_upper)==as.character(ta.t0)]
        ta_lb <- ci_lower[names(ci_lower)==as.character(ta.t0)]
      } else {
        ta_ub <- betahat[which(timeVec==ta.t0)]+ta.s*sqrt(diag(covhat)[which(timeVec==ta.t0)])
        ta_lb <- betahat[which(timeVec==ta.t0)]-ta.s*sqrt(diag(covhat)[which(timeVec==ta.t0)])
      }

      honest_ub_splinefun <- function(x) ta_ub
      honest_lb_splinefun <- function(x) ta_lb

      honest_ub_splinefun <- Vectorize(honest_ub_splinefun)
      honest_lb_splinefun <- Vectorize(honest_lb_splinefun)

    }


    if(!is.null(frm.mbar) & is.null(ftr.m) & is.null(frmtr.mbar)) {

      if(is.null(ta.s)) {
        ta_ub <- ci_upper[names(ci_upper)==as.character(ta.t0)]
        ta_lb <- ci_lower[names(ci_lower)==as.character(ta.t0)]
      } else {
        ta_ub <- betahat[which(timeVec==ta.t0)]+ta.s*sqrt(diag(covhat)[which(timeVec==ta.t0)])
        ta_lb <- betahat[which(timeVec==ta.t0)]-ta.s*sqrt(diag(covhat)[which(timeVec==ta.t0)])
      }

      # frm bounds
      ts             <- seq(start, ta.t0, length.out=50*sum(timeVec<=ta.t0))
      betas_deriv    <- pracma::fderiv(betahat_splinefun, ts, n=1, method="backward")
      mean_abs_deriv <- mean(abs(betas_deriv))

      frm_ub_ta_ub <- function (x) {ifelse(x>=ta.t0, frm.mbar*mean_abs_deriv*(x-ta.t0)+ta_ub, NA)}
      frm_lb_ta_ub <- function (x) {ifelse(x>=ta.t0, -frm.mbar*mean_abs_deriv*(x-ta.t0)+ta_ub, NA)}
      frm_ub_ta_lb <- function (x) {ifelse(x>=ta.t0, frm.mbar*mean_abs_deriv*(x-ta.t0)+ta_lb, NA)}
      frm_lb_ta_lb <- function (x) {ifelse(x>=ta.t0, -frm.mbar*mean_abs_deriv*(x-ta.t0)+ta_lb, NA)}

      x_vals <- seq(start, end, length.out = 5*length(timeVec))

      if (ta.t0!=t0) {
        honest_ub_vals <- apply(cbind(ta_ub, frm_ub_ta_ub(x_vals), frm_ub_ta_lb(x_vals)), 1,  function(row) max(row,na.rm=TRUE))
        honest_lb_vals <- apply(cbind(ta_lb, frm_lb_ta_ub(x_vals), frm_lb_ta_lb(x_vals)), 1,  function(row) min(row,na.rm=TRUE))
      } else {
        honest_ub_vals <- apply(cbind(frm_ub_ta_ub(x_vals), frm_ub_ta_lb(x_vals)), 1,  function(row) if (all(is.na(row))) 0 else max(row, na.rm=TRUE))
        honest_lb_vals <- apply(cbind(frm_lb_ta_ub(x_vals), frm_lb_ta_lb(x_vals)), 1,  function(row) if (all(is.na(row))) 0 else min(row,na.rm=TRUE))
      }
      honest_ub_splinefun <- splinefun(x=x_vals, y=honest_ub_vals, method="natural")
      honest_lb_splinefun <- splinefun(x=x_vals, y=honest_lb_vals, method="natural")
      honest_ub_splinefun <- Vectorize(honest_ub_splinefun)
      honest_lb_splinefun <- Vectorize(honest_lb_splinefun)
    }


    if(is.null(frm.mbar) & !is.null(ftr.m) & is.null(frmtr.mbar)) {

      if(is.null(ta.s)) {
        ta_ub <- ci_upper[names(ci_upper)==as.character(ta.t0)]
        ta_lb <- ci_lower[names(ci_lower)==as.character(ta.t0)]
      } else {
        ta_ub <- betahat[which(timeVec==ta.t0)]+ta.s*sqrt(diag(covhat)[which(timeVec==ta.t0)])
        ta_lb <- betahat[which(timeVec==ta.t0)]-ta.s*sqrt(diag(covhat)[which(timeVec==ta.t0)])
      }

      # ftr bounds
      ts          <- seq(start, ta.t0, length.out=50*sum(timeVec<=ta.t0))
      betas_deriv <- pracma::fderiv(betahat_splinefun, ts, n=1, method="backward")
      slope       <- mean(betas_deriv)

      ftr_ub_ta_ub <- function (x) {ifelse(x>=ta.t0, (ftr.m*abs(slope)+slope)*(x-ta.t0)+ta_ub, NA)}
      ftr_lb_ta_ub <- function (x) {ifelse(x>=ta.t0, (-ftr.m*abs(slope)+slope)*(x-ta.t0)+ta_ub, NA)}
      ftr_ub_ta_lb <- function (x) {ifelse(x>=ta.t0, (ftr.m*abs(slope)+slope)*(x-ta.t0)+ta_lb, NA)}
      ftr_lb_ta_lb <- function (x) {ifelse(x>=ta.t0, (-ftr.m*abs(slope)+slope)*(x-ta.t0)+ta_lb, NA)}

      x_vals <- seq(start, end, length.out=5*length(timeVec))

      if (ta.t0!=t0) {
        honest_ub_vals <- apply(cbind(ta_ub, ftr_ub_ta_ub(x_vals), ftr_ub_ta_lb(x_vals)), 1,  function(row) max(row,na.rm=TRUE))
        honest_lb_vals <- apply(cbind(ta_lb, ftr_lb_ta_ub(x_vals), ftr_lb_ta_lb(x_vals)), 1,  function(row) min(row,na.rm=TRUE))
      } else {
        honest_ub_vals <- apply(cbind(ftr_ub_ta_ub(x_vals), ftr_ub_ta_lb(x_vals)), 1,  function(row) if (all(is.na(row))) 0 else max(row,na.rm=TRUE))
        honest_lb_vals <- apply(cbind(ftr_lb_ta_ub(x_vals), ftr_lb_ta_lb(x_vals)), 1,  function(row) if (all(is.na(row))) 0 else min(row,na.rm=TRUE))
      }
      honest_ub_splinefun <- splinefun(x=x_vals, y=honest_ub_vals, method="natural")
      honest_lb_splinefun <- splinefun(x=x_vals, y=honest_lb_vals, method="natural")
      honest_ub_splinefun <- Vectorize(honest_ub_splinefun)
      honest_lb_splinefun <- Vectorize(honest_lb_splinefun)
    }

    if(is.null(frm.mbar) & is.null(ftr.m) & !is.null(frmtr.mbar)) {

      if(is.null(ta.s)) {
        ta_ub <- ci_upper[names(ci_upper)==as.character(ta.t0)]
        ta_lb <- ci_lower[names(ci_lower)==as.character(ta.t0)]
      } else {
        ta_ub <- betahat[which(timeVec==ta.t0)]+ta.s*sqrt(diag(covhat)[which(timeVec==ta.t0)])
        ta_lb <- betahat[which(timeVec==ta.t0)]-ta.s*sqrt(diag(covhat)[which(timeVec==ta.t0)])
      }

      # frmtr bounds
      ts             <- seq(start, ta.t0, length.out=50*sum(timeVec<=ta.t0))
      betas_deriv    <- pracma::fderiv(betahat_splinefun, ts, n=1, method="backward")
      mean_abs_deriv <- mean(abs(betas_deriv))
      slope          <- mean(betas_deriv)

      ftr_ub_ta_ub <- function (x) {ifelse(x>=ta.t0, (frmtr.mbar*mean_abs_deriv+slope)*(x-ta.t0)+ta_ub, NA)}
      ftr_lb_ta_ub <- function (x) {ifelse(x>=ta.t0, (-frmtr.mbar*mean_abs_deriv+slope)*(x-ta.t0)+ta_ub, NA)}
      ftr_ub_ta_lb <- function (x) {ifelse(x>=ta.t0, (frmtr.mbar*mean_abs_deriv+slope)*(x-ta.t0)+ta_lb, NA)}
      ftr_lb_ta_lb <- function (x) {ifelse(x>=ta.t0, (-frmtr.mbar*mean_abs_deriv+slope)*(x-ta.t0)+ta_lb, NA)}

      x_vals <- seq(start, end, length.out=5*length(timeVec))

      if (ta.t0!=t0) {
        honest_ub_vals <- apply(cbind(ta_ub, ftr_ub_ta_ub(x_vals), ftr_ub_ta_lb(x_vals)), 1,  function(row) max(row,na.rm=TRUE))
        honest_lb_vals <- apply(cbind(ta_lb, ftr_lb_ta_ub(x_vals), ftr_lb_ta_lb(x_vals)), 1,  function(row) min(row,na.rm=TRUE))
      } else {
        honest_ub_vals <- apply(cbind(ftr_ub_ta_ub(x_vals), ftr_ub_ta_lb(x_vals)), 1,  function(row) if (all(is.na(row))) 0 else max(row,na.rm=TRUE))
        honest_lb_vals <- apply(cbind(ftr_lb_ta_ub(x_vals), ftr_lb_ta_lb(x_vals)), 1,  function(row) if (all(is.na(row))) 0 else min(row,na.rm=TRUE))
      }
      honest_ub_splinefun <- splinefun(x=x_vals, y=honest_ub_vals, method="natural")
      honest_lb_splinefun <- splinefun(x=x_vals, y=honest_lb_vals, method="natural")
      honest_ub_splinefun <- Vectorize(honest_ub_splinefun)
      honest_lb_splinefun <- Vectorize(honest_lb_splinefun)
    }

    # find the time spans of statistical significance
    n.int <-200

    roots_UB <- vector("numeric", length = n.int)
    for (i in 1:n.int) {
      fun         <- function(x) (scb_ub_splinefun(x)-honest_ub_splinefun(x))*(scb_lb_splinefun(x)-honest_ub_splinefun(x))
      interval    <- c(start+(i-1)*(end-start)/n.int, start+(i)*(end-start)/n.int)
      root        <- try(uniroot(fun, interval = interval)$root, silent = TRUE)
      roots_UB[i] <- if(inherits(root, "try-error")) {NA} else root
    }
    roots_UB <- sort(unique(roots_UB[!is.na(roots_UB)]))

    roots_LB <- vector("numeric", length = n.int)
    for (i in 1:n.int) {
      fun         <- function(x) (scb_ub_splinefun(x)-honest_lb_splinefun(x))*(scb_lb_splinefun(x)-honest_lb_splinefun(x))
      interval    <- c(start+(i-1)*(end-start)/n.int, start+(i)*(end-start)/n.int)
      root        <- try(uniroot(fun, interval = interval)$root, silent = TRUE)
      roots_LB[i] <- if(inherits(root, "try-error")) {NA} else root
    }
    roots_LB <-  sort(unique(roots_LB[!is.na(roots_LB)]))

    roots_UB_vec <- sort(unique(c(start,roots_UB,end)))
    roots_LB_vec <- sort(unique(c(start,roots_LB,end)))
    roots_vec    <- sort(unique(c(start,end,roots_UB, roots_LB)))

    stat_sig_UB_vec <- rep(NA, length(roots_UB_vec)-1)
    stat_sig_LB_vec <- rep(NA, length(roots_LB_vec)-1)

    for (i in 1:(length(roots_UB_vec)-1)) {
      test_point <- (roots_UB_vec[i]+roots_UB_vec[i+1])/2
      scb_ub_sign <- sign(scb_ub_splinefun( test_point)-honest_ub_splinefun(test_point))
      scb_lb_sign <- sign(scb_lb_splinefun( test_point)-honest_ub_splinefun(test_point))
      stat_sig_UB_vec[i] <- (scb_ub_sign+scb_lb_sign)/2
    }

    for (i in 1:(length(roots_LB_vec)-1)) {
      test_point <- (roots_LB_vec[i]+roots_LB_vec[i+1])/2
      scb_ub_sign <- sign(scb_ub_splinefun( test_point)-honest_lb_splinefun(test_point))
      scb_lb_sign <- sign(scb_lb_splinefun( test_point)-honest_lb_splinefun(test_point))
      stat_sig_LB_vec[i] <- (scb_ub_sign+scb_lb_sign)/2
    }

    interval_fun <- function(x, a, b) {
      if (length(b) != length(a) - 1) {
        stop("Vector 'b' must be one element shorter than vector 'a'.")
      }
      for (i in 1:(length(a) - 1)) {
        if (x >= a[i] && x <= a[i + 1]) {
          return(b[i])
        }
      }
      return(NA)
    }

    fun_stat_sig_UB_vec <- function(x) interval_fun(x, roots_UB_vec, stat_sig_UB_vec)
    fun_stat_sig_LB_vec <- function(x) interval_fun(x, roots_LB_vec, stat_sig_LB_vec)
    fun_stat_sig_UB_vec <- Vectorize(fun_stat_sig_UB_vec)
    fun_stat_sig_LB_vec <- Vectorize(fun_stat_sig_LB_vec)

    fun_stat_sig_vec <- function(x) {if(fun_stat_sig_UB_vec(x) == fun_stat_sig_LB_vec(x)) fun_stat_sig_UB_vec(x) else 0}
    fun_stat_sig_vec <- Vectorize(fun_stat_sig_vec)

    # draw plots
    x_vals  <- seq(start, end, length.out=5*length(timeVec))
    y_range <- range(c(betahat_splinefun(x_vals), scb_ub_splinefun(x_vals), scb_lb_splinefun(x_vals)), na.rm=TRUE)

    plot(timeVec, betahat,  type = "p", ylim =y_range, xlab = "Event Time", ylab = "", pch = 16, col = "blue", cex.axis=1.4, cex.lab=1.4, ...)
    title(ylab=expression(hat(beta)), mgp=c(2.1,1,0), cex.lab=1.4)
    grid()

    if(isTRUE(ci)) arrows(timeVec, ci_lower, timeVec, ci_upper, angle = 90, code = 3, length = 0.025, col = "blue4")
    curve(betahat_splinefun, add=TRUE, lty=1)
    curve(scb_ub_splinefun, add=TRUE, lty=1)
    curve(scb_lb_splinefun, add=TRUE, lty=1)

    ShadeBetween <- function(x1, x2, f1, f2, ...) {polygon(c(x1, rev(x2)), c(f1, rev(f2)), ...)}

    if(isTRUE(post.trt)) {
      segments(x0=t0, y0=0, x1=end, y1=0, lty=3, lwd=4, col="red")
      ShadeBetween(timeVec[timeVec>=t0], timeVec[timeVec>=t0], honest_ub_splinefun(timeVec[timeVec>=t0]), honest_lb_splinefun(timeVec[timeVec>=t0]), col=rgb(1,0,0,alpha=0.3), border=NA)

      if(!(t0 %in%  roots_vec)) {roots_vec <- sort(c(roots_vec, t0))}

      for (i in 1:(length(roots_vec)-1) ) {
        if ( fun_stat_sig_vec((roots_vec[i]+roots_vec[i+1])/2)==1 & roots_vec[i]>=t0 ) {
          rect(roots_vec[i], y_range[1], roots_vec[i+1], y_range[2], col=rgb(0.5, 0.5, 0.5, alpha=0.7), border=NA)
        } else {
          if (fun_stat_sig_vec((roots_vec[i]+roots_vec[i+1])/2)==-1 & roots_vec[i]>=t0 ) {
            rect(roots_vec[i], y_range[1], roots_vec[i+1], y_range[2], col=rgb(0.5, 0.5, 0.5, alpha=0.4), border=NA)
          }
        }
      }
    } else {
      segments(x0=start, y0=0, x1=end, y1=0, lty=3, lwd=4, col="red")
      ShadeBetween(timeVec, timeVec, honest_ub_splinefun(timeVec), honest_lb_splinefun(timeVec), col=rgb(1,0,0,alpha=0.3), border=NA)

      for (i in 1:(length(roots_vec)-1) ) {
        if ( fun_stat_sig_vec((roots_vec[i]+roots_vec[i+1])/2)==1  ) {
          rect(roots_vec[i], y_range[1], roots_vec[i+1], y_range[2], col=rgb(0.5, 0.5, 0.5, alpha=0.7), border=NA)
        } else {
          if (fun_stat_sig_vec((roots_vec[i]+roots_vec[i+1])/2)==-1) {
            rect(roots_vec[i], y_range[1], roots_vec[i+1], y_range[2], col=rgb(0.5, 0.5, 0.5, alpha=0.4), border=NA)
          }
        }
      }
    }

    # write a legend for the plot
    if(!is.null(pos.legend)) {

      xlim <- par("usr")[1:2]
      ylim <- par("usr")[3:4]

      if (ta.t0!=t0 & (!is.null(frm.mbar) | !is.null(ftr.m) | !is.null(frmtr.mbar))) {upper_text <- "Reference values with anticipation and PT violations"; lower_text <- "Reference value without anticipation or PT violations          "}
      if (ta.t0==t0 & (!is.null(frm.mbar) | !is.null(ftr.m) | !is.null(frmtr.mbar))) {upper_text <- "Reference values with PT violations"; lower_text <- "Reference value without PT violations          "}
      if (ta.t0!=t0 & is.null(frm.mbar) & is.null(ftr.m) & is.null(frmtr.mbar)) {upper_text <- "Reference values with anticipation"; lower_text <- "Reference value without anticipation          "}

      text_width_upper   <- strwidth(upper_text) # calculate the width of the text dynamically
      text_width_lower   <- strwidth(lower_text)
      fixed_rect_b_width <- strwidth("Reference") # define a longer fixed width for the red rectangle (B) and red dashed line
      fixed_line_width   <- strwidth("Reference")
      extra_padding      <- strwidth("   ", cex=1) # add extra padding for the white rectangle to avoid text touching borders
      rect_width         <- max(text_width_upper, text_width_lower)+extra_padding+fixed_rect_b_width  # calculate the width of the white rectangle (A) based on the longest text
      rect_height        <- diff(ylim)*0.1 # define the height for the white rectangle (A)
      rect_x1            <- mean(xlim)-rect_width/2 # define positions for the white rectangle A (centered)
      rect_x2            <- mean(xlim)+rect_width/2

      if(pos.legend=="top")    {rect_y2 <- ylim[2]-rect_height*0.5; rect_y1 <- rect_y2-rect_height}
      if(pos.legend=="bottom") {rect_y2 <- ylim[1]+rect_height*0.5; rect_y1 <- rect_y2+rect_height}

      rect(rect_x1, rect_y1, rect_x2, rect_y2, col = "white", border = "black") # draw white rectangular A with black border

      left_padding <- rect_width*0.05 # define padding to prevent touching the left border
      rect_b_x1    <- rect_x1+left_padding  # position and draw the fixed-width red rectangle B (aligned with upper text)
      rect_b_x2    <- rect_b_x1+fixed_rect_b_width

      if(pos.legend=="top")    {rect_b_y2 <- rect_y2-rect_height*0.15; rect_b_y1 <- rect_b_y2-rect_height*0.35}
      if(pos.legend=="bottom") {rect_b_y2 <- rect_y1-rect_height*0.5; rect_b_y1 <- rect_b_y2+rect_height*0.35}

      rect(rect_b_x1, rect_b_y1, rect_b_x2, rect_b_y2, col = rgb(1, 0, 0, alpha = 0.3), border = NA)  # draw red rectangle B (transparent, no border)
      text(rect_b_x2 + left_padding, (rect_b_y1 + rect_b_y2) / 2, upper_text, cex = 1,  adj = 0)  # add upper text next to red rectangle B

      line_x1 <- rect_b_x1  # define position for the red dashed line (aligned with lower text) # same as red rectangle start
      line_x2 <- rect_b_x2  # same fixed width as red rectangle

      if(pos.legend=="top")    {line_y <- rect_y1+rect_height*0.3}
      if(pos.legend=="bottom") {line_y <- rect_y2+rect_height*0.3}

      segments(line_x1, line_y, line_x2, line_y, col = "red", lwd = 4, lty = 3)  # draw the red dashed line (fixed width)
      text(line_x2 + left_padding, line_y, lower_text, cex = 1, adj = 0)  # add lower text next to the red dashed line
    }
  }
  options(warn=1)
}
