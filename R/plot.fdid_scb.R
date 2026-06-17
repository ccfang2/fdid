#' Custom plot method for objects of class \code{"fdid_scb"}
#' @description This function is an S3 method for \link[base]{plot}, specifically designed for objects of class \code{"fdid_scb"}. It plots the simultaneous confidence bands together with honest reference band for performing honest causal inference.
#'
#' @param object an object of class \code{"fdid_scb"}. The object to be plotted.
#' @param ta.ts a numeric value that indicates the last time point in the pre-anticipation periods. It should be smaller than the reference time period.
#' @param ta.s a numeric vector of control parameters with length two for deriving the honest reference band under the violation of no-anticipation assumption. If \code{ta.ts} is NULL, the value of \code{ta.s} is ignored. If \code{ta.ts} is not NULL, the value of \code{ta.s} cannot be NULL.
#' @param ftr.m a numeric vector of control parameters with length two for deriving the honest reference band under violation of the parallel trends assumption. The default is NULL.
#' @param ref.band.pre a logical value. If TRUE, the reference band for pre-anticipation period is also plotted.
#' @param note.pre a logical value. If TRUE, the note for pre-anticipation period is given on top of plot. If no honest reference band is defined, there is no need to perform validation in the pre-anticipation period, and \code{note.pre} is FALSE.
#' @param note.post a logical value. If TRUE, the note for post-treatment period is given on top of plot.
#' @param ci.pre a logical value. If TRUE, the point-wise confidence intervals for pre-treatment period are plotted.
#' @param ci.post a logical value. If TRUE, the point-wise confidence intervals for post-treatment period are plotted.
#' @param pos.legend a character value of "top" or "bottom" that indicates the position of legend. If NULL, the legend is not printed.
#' @param scale.legend a positive number that defines the size of legend. If \code{pos.legend} is NULL, the value of \code{scale.legend} is ignored.
#' @param verbose a logical value. If TRUE, the uniformly significant time span(s) under your specification will be directly printed out on the console.
#' @param ... Additional arguments to be passed to \link[base]{plot}.
#'
#' @return The function returns a plot with simultaneous confidence bands for event study coefficients in a functional framework, together with the honest reference band for honest inference, if properly defined.
#'
#' @import pracma
#' @export
#'
#' @references Fang, C. and Liebl, D. (2026). Making Event Study Plots Honest: A Functional Data Approach to Causal Inference. \href{https://arxiv.org/abs/2512.06804}{arXiv:2512.06804}.
#' @seealso \link{fdid_scb}
#'
#' @examples
#' data(LWdata)
#' fdid_scb_est <- fdid_scb(beta=LWdata$beta, cov=LWdata$cov, t0=LWdata$t0)
#' cat("The reference time is ", LWdata$t0, ". If not NULL, the input 'ta.ts' in function 'plot' should be smaller than this value.", sep="")
#'
#' ## simultaneous inference
#' par(cex.axis = 1.4, cex.lab = 1.4, cex.main = 1.4)
#' plot(fdid_scb_est, scale.legend=1.4)
#' title(ylab="Effects of Duty-to-Bargain Laws")
#'
#' ## honest inference under treatment anticipation
#' plot(fdid_scb_est, ta.ts=-3, ta.s=c(0.2,0.2), scale.legend=1.4)
#'
#' ## honest inference under violation of parallel trends assumption
#' plot(fdid_scb_est, ftr.m=c(0.2,0.2), scale.legend=1.4)
plot.fdid_scb <- function(object,
                          ta.ts=NULL,
                          ta.s=NULL,
                          ftr.m=NULL,
                          ref.band.pre=TRUE,
                          note.pre=TRUE,
                          note.post=TRUE,
                          ci.pre=FALSE,
                          ci.post=FALSE,
                          pos.legend="top",
                          scale.legend=1,
                          verbose=TRUE,
                          ...) {

  options(warn=-1)

  # check conditions
  if (!base::inherits(object,"fdid_scb")) stop("The input 'object' should be an output of function 'fdid_scb'.")
  if (!is.null(ta.ts) && (!is.numeric(ta.ts) || length(ta.ts) != 1)) stop("The input 'ta.ts' should be either NULL or a numeric scalar.")
  if (!is.null(ta.ts) && !(ta.ts %in% object$scb$event_t[which(object$scb$event_t < object$data$t0)])) stop("If not NULL, the input 'ta.ts' should be among the pre-treatment event time in 'object', excluding the reference time point.")
  if (!is.null(ta.ts) && ta.ts == object$scb$event_t[1] && !all(sapply(list( ftr.m), is.null))) stop("If 'ta.ts' is defined to be the first event time, 'ftr.m' must be NULL, because there is no available data for computing pre-trend differences.")
  if (!is.null(ta.s) && any(is.na(ta.s))) stop("The input 'ta.s' cannot contain NA.")
  if (!is.null(ta.s) && (!is.numeric(ta.s) || length(ta.s) != 2 || any(ta.s <0) )) stop("The input 'ta.s' should be either NULL or a numeric non-negative vector of length two.")
  if (!is.null(ta.ts) && is.null(ta.s)) stop("If the input 'ta.ts' is not NULL, 'ta.s' should not be NULL.")
  if (!is.null(ftr.m) && any(is.na(ftr.m))) stop("The input 'ftr.m' cannot contain NA.")
  if (!is.null(ftr.m) && (!is.numeric(ftr.m) || length(ftr.m) != 2 || any(ftr.m<0) )) stop("The input 'ftr.m' should be either NULL or a numeric non-negative vector of length two.")
  if (!is.logical(ref.band.pre)) stop("The input 'ref.band.pre' should be logical.")
  if (!is.logical(note.pre)) stop("The input 'note.pre' should be logical.")
  if (!is.logical(note.post)) stop("The input 'note.post' should be logical.")
  if (!is.logical(ci.pre)) stop("The input 'ci.pre' should be logical.")
  if (!is.logical(ci.post)) stop("The input 'ci.post' should be logical.")
  if (!is.null(pos.legend) && !pos.legend %in% c("top", "bottom")) stop("The input 'pos.legend' must be 'top', 'bottom', or NULL.")
  if (!is.numeric(scale.legend) || length(scale.legend) != 1 || scale.legend <= 0 || !is.finite(scale.legend)) stop("The input 'scale.legend' must be a positive number.")
  if (!is.logical(verbose)) stop("The input 'verbose' should be logical.")

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

  ci.alpha                  <- object$data$ci.alpha
  scb.pre.alpha             <- object$data$scb.pre.alpha
  scb.post.alpha            <- object$data$scb.post.alpha

  # regular statistical inference
  if (is.null(ta.ts) & is.null(ftr.m)) {

    t0_m1 <- timeVec[which(timeVec==t0)-1]
    t0_p1 <- timeVec[which(timeVec==t0)+1]

    # find the time spans of statistical significance
    n.int <- 200
    roots <- vector("numeric", length = n.int)
    for (i in 1:n.int) {
      fun      <- function(x) scb_ub_splinefun(x)*scb_lb_splinefun(x)
      interval <- c(start+(i-1)*(end-start)/n.int, start+(i)*(end-start)/n.int)
      root     <- try(uniroot(fun, interval = interval)$root, silent = TRUE)
      roots[i] <- if(inherits(root, "try-error")) {NA} else root
    }

    roots <- sort(unique(roots[!is.na(roots)]))

    roots_vec    <- sort(unique(c(start,roots,end)))
    stat_sig_vec <- rep(NA, length(roots_vec)-1)
    for (i in 1:(length(roots_vec)-1)) {
      test_point      <- (roots_vec[i]+roots_vec[i+1])/2
      scb_ub_sign     <- sign(scb_ub_splinefun( test_point))
      scb_lb_sign     <- sign(scb_lb_splinefun( test_point))
      stat_sig_vec[i] <- (scb_ub_sign+scb_lb_sign)/2
    }

    x_vals  <- seq(start, end, length.out = 7*length(timeVec))

    # draw plot
    y_range     <- range(c(betahat_splinefun(x_vals), scb_ub_splinefun(x_vals), scb_lb_splinefun(x_vals)), na.rm=TRUE)
    y_range     <- c(y_range[1] - diff(y_range)*0.12, y_range[2] + diff(y_range)*0.08)
    y_range_fun <- function(x) y_range[1]+(y_range[2]-y_range[1])*(x-start)/(end-start)
    curve(y_range_fun, from=start, to=end, lty=1,  xlab = "Event Time", ylab = "", col="white",n=500,...) #ylim =y_range,

    ShadeBetween <- function(x1, x2, f1, f2, ...) {polygon(c(x1, rev(x2)), c(f1, rev(f2)), ...)}
    band_grid    <- seq(start, end, length.out = 400)

    ub <- scb_ub_splinefun(seq(start, end, length.out = 400))
    lb <- scb_lb_splinefun(seq(start, end, length.out = 400))

    polygon(
      x = c(seq(t0, end, length.out = 200), rev(seq(t0, end, length.out = 200))),
      y = c(scb_ub_splinefun(seq(t0, end, length.out = 200)), rev(scb_lb_splinefun(seq(t0, end, length.out = 200)))),
      col =  adjustcolor("#48C9B0", alpha.f = 1),
      border = NA
      )

    polygon(
      x = c(seq(start, t0, length.out = 200), rev(seq(start, t0, length.out = 200))),
      y = c(scb_ub_splinefun(seq(start, t0, length.out = 200)), rev(scb_lb_splinefun(seq(start, t0, length.out = 200)))),
      col = adjustcolor("#48C9B0", alpha.f = 0.7),
      border = NA
      )

    curve(betahat_splinefun, from=start, to=end, lty=1, add=TRUE, xlab = "Event Time", ylab = "", col="blue", n=500, ...) #ylim =y_range,
    points(timeVec[which(timeVec<=t0)], betahat[which(timeVec<=t0)], pch = 16, col = rgb(0,0,1,alpha=0.4), cex=1)
    points(timeVec[which(timeVec>t0)], betahat[which(timeVec>t0)], pch = 16, col = "blue", cex=1)

    if(isTRUE(ci.pre)) {arrows(timeVec[which(timeVec<=t0)], ci_lower[which(timeVec<=t0)], timeVec[which(timeVec<=t0)], ci_upper[which(timeVec<=t0)], angle = 90, code = 3, length = 0.025, col = rgb(0,0,1,alpha=0.4));}
    if(isTRUE(ci.post)) {arrows(timeVec[which(timeVec>t0)], ci_lower[which(timeVec>t0)], timeVec[which(timeVec>t0)], ci_upper[which(timeVec>t0)], angle = 90, code = 3, length = 0.025, col = "blue4");}

    abline(v=(t0+timeVec[which(timeVec==t0)+1])/2, lty=3, col="blue4", lwd=4)

    # significant intervals
    if(!isTRUE(ref.band.pre)) {
      # segments(x0=t0, y0=0, x1=end, y1=0, lty=3, lwd=4, col="red")
      segments(x0=t0_p1, y0=0, x1=end, y1=0, lty=3, lwd=4, col="red")
      roots_vec_temp    <- sort(unique(c(roots_vec, timeVec[which(timeVec==t0_p1)])))
      intervals         <- data.frame(left  = head(roots_vec_temp, -1), right = tail(roots_vec_temp, -1))
      stat_sig_vec_temp <- rep(NA, nrow(intervals))

      for (i in 1:nrow(intervals)) {
        l <- intervals$left[i]
        r <- intervals$right[i]
          idx <- which(roots_vec <= l)[length(which(roots_vec  <= l))]
          stat_sig_vec_temp[i] <- stat_sig_vec[idx]
          }

      stat_sig_vec <- stat_sig_vec_temp
      roots_vec    <- roots_vec_temp

      for (i in 1:(length(roots_vec)-1)) {
        if (stat_sig_vec[i]==1 & roots_vec[i]>=t0_p1 ) {
          if(isTRUE(verbose)) cat("Uniformly and Positively Significant Time Span: [", roots_vec[i], ",", roots_vec[i+1], "]\n")
        } else {
          if (stat_sig_vec[i]==-1 & roots_vec[i]>=t0_p1) {
            if(isTRUE(verbose)) cat("Uniformly and Negatively Significant Time Span: [", roots_vec[i], ",", roots_vec[i+1], "]\n")
          }
        }
      }
    } else {
      segments(x0=start, y0=0, x1=end, y1=0, lty=3, lwd=4, col="red")
      roots_vec_temp    <- sort(unique(c(roots_vec, timeVec[which(timeVec==t0_p1)])))
      intervals         <- data.frame(left  = head(roots_vec_temp, -1), right = tail(roots_vec_temp, -1))
      stat_sig_vec_temp <- rep(NA, nrow(intervals))

      for (i in 1:nrow(intervals)) {
        l <- intervals$left[i]
        r <- intervals$right[i]
          idx <- which(roots_vec <= l)[length(which(roots_vec  <= l))]
          stat_sig_vec_temp[i] <- stat_sig_vec[idx]
          }

      stat_sig_vec <- stat_sig_vec_temp
      roots_vec    <- roots_vec_temp

      for (i in 1:(length(roots_vec)-1)) {
        if (stat_sig_vec[i]==1 & roots_vec[i]>=t0_p1 ) {
          if (isTRUE(verbose)) cat("Uniformly and Positively Significant Time Span: [", roots_vec[i], ",", roots_vec[i+1], "]\n")
        } else {
          if (stat_sig_vec[i]==-1 & roots_vec[i]>=t0_p1 ) {
            if (isTRUE(verbose)) cat("Uniformly and Negatively Significant Time Span: [", roots_vec[i], ",", roots_vec[i+1], "]\n")
          }
        }
      }
    }

    # write a legend for the plot
    if(!is.null(pos.legend)) {

      xlim <- par("usr")[1:2]
      ylim <- par("usr")[3:4]

      text              <- "Classical Reference Line          "
      text_width        <- strwidth(text, cex = scale.legend)
      fixed_line_width  <- strwidth("Classical", cex=scale.legend)
      extra_padding     <- strwidth("   ", cex=scale.legend)
      rect_width        <- text_width+extra_padding+fixed_line_width
      rect_height       <- diff(ylim)*0.08*scale.legend
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

      segments(line_x1, line_y, line_x2, line_y, col = "red", lwd = 4*scale.legend, lty = 3)
      text(line_x2 + left_padding, line_y, text, cex = scale.legend, adj = 0)
    }

    if (isTRUE(note.post)) {

      x_left   <- par("usr")[1]
      x_right  <- par("usr")[2]
      y_top    <- par("usr")[4]
      y_bottom <- par("usr")[3]

      box_height <- 0.1 * (y_top - y_bottom)  # height = 10% of plot range
      x_mid      <- timeVec[which(timeVec==t0)+1]      # border between boxes

      par(xpd = TRUE) # Turn off clipping so we can draw *outside* the plot
      segments(x_mid,  y_top,  x_mid,  y_top + box_height, col = "black")

      text((x_mid + x_right)/2, y_top + 0.5*box_height, "Testing", cex = scale.legend)
      par(xpd = FALSE) # Turn clipping back on
    }
  }

  # honest inference

  # ta bounds
    if( !is.null(ta.ts) & is.null(ftr.m)) {

      t0_p1 <- timeVec[which(timeVec==t0)+1]

      honest_ub_splinefun <- function(x) mean(betahat[which(timeVec<=ta.ts)]) + ta.s[1]
      honest_lb_splinefun <- function(x) mean(betahat[which(timeVec<=ta.ts)]) - ta.s[2]

      honest_ub_splinefun <- Vectorize(honest_ub_splinefun)
      honest_lb_splinefun <- Vectorize(honest_lb_splinefun)
    }

  # ftr bounds
  if( is.null(ta.ts) & !is.null(ftr.m)) {

    t0_m1 <- timeVec[which(timeVec==t0)-1]
    t0_p1 <- timeVec[which(timeVec==t0)+1]

    betahat_slopes <- sapply(1:which(timeVec==t0_m1), function(index) (betahat[index]-betahat[which(timeVec==t0)])/(timeVec[index]-t0))

    honest_ub_splinefun <- function(x) mean(betahat_slopes)*(t0_p1-t0)+ftr.m[1]+mean(betahat_slopes)*(x-t0_p1)
    honest_lb_splinefun <- function(x) mean(betahat_slopes)*(t0_p1-t0) - ftr.m[2] + mean(betahat_slopes)*(x-t0_p1)

    honest_ub_splinefun <- Vectorize(honest_ub_splinefun)
    honest_lb_splinefun <- Vectorize(honest_lb_splinefun)
  }

  # ta + ftr bounds
    if(!is.null(ta.ts) & !is.null(ftr.m)) {

      ta_ub <- mean(betahat[which(timeVec<=ta.ts)]) + ta.s[1]
      ta_lb <- mean(betahat[which(timeVec<=ta.ts)]) - ta.s[2]

      ta.ts_m1 <- timeVec[which(timeVec==ta.ts)-1]
      ta.ts_p1 <- timeVec[which(timeVec==ta.ts)+1]

      betahat_slopes <- sapply(1:which(timeVec==ta.ts_m1), function(index) (betahat[index]-betahat[which(timeVec==ta.ts)])/(timeVec[index]-ta.ts))

      ftr_ub_ta_ub <- function (x) ftr.m[1]+ mean(betahat_slopes)*(x-ta.ts)+ta_ub
      ftr_lb_ta_ub <- function (x) -ftr.m[2]+ mean(betahat_slopes)*(x-ta.ts)+ta_ub
      ftr_ub_ta_lb <- function (x) ftr.m[1]+ mean(betahat_slopes)*(x-ta.ts)+ta_lb
      ftr_lb_ta_lb <- function (x) -ftr.m[2]+ mean(betahat_slopes)*(x-ta.ts)+ta_lb

      x_vals <- seq(start, end, length.out=7*length(timeVec))

      honest_ub_vals <- apply(cbind(ta_ub, ftr_ub_ta_ub(x_vals), ftr_ub_ta_lb(x_vals)), 1,  function(row) max(row,na.rm=TRUE))
      honest_lb_vals <- apply(cbind(ta_lb, ftr_lb_ta_ub(x_vals), ftr_lb_ta_lb(x_vals)), 1,  function(row) min(row,na.rm=TRUE))
      honest_ub_splinefun <- splinefun(x=x_vals, y=honest_ub_vals, method="natural")
      honest_lb_splinefun <- splinefun(x=x_vals, y=honest_lb_vals, method="natural")
      honest_ub_splinefun <- Vectorize(honest_ub_splinefun)
      honest_lb_splinefun <- Vectorize(honest_lb_splinefun)
    }

  if(!(is.null(ta.ts) & is.null(ftr.m))) {

    # find the time spans of statistical significance
    n.int <-200

    t0_p1    <- timeVec[which(timeVec==t0)+1]
    roots_UB <- vector("numeric", length = n.int)

    for (i in 1:n.int) {
      fun         <- function(x) (scb_ub_splinefun(x)-honest_ub_splinefun(x))*(scb_lb_splinefun(x)-honest_ub_splinefun(x))
      interval    <- c(t0_p1+(i-1)*(end-t0_p1)/n.int, t0_p1+(i)*(end-t0_p1)/n.int)
      root        <- try(uniroot(fun, interval = interval)$root, silent = TRUE)
      roots_UB[i] <- if(inherits(root, "try-error")) {NA} else root
    }
    roots_UB <- sort(unique(roots_UB[!is.na(roots_UB)]))

    roots_LB <- vector("numeric", length = n.int)
    for (i in 1:n.int) {
      fun         <- function(x) (scb_ub_splinefun(x)-honest_lb_splinefun(x))*(scb_lb_splinefun(x)-honest_lb_splinefun(x))
      interval    <- c(t0_p1+(i-1)*(end- t0_p1)/n.int,  t0_p1+(i)*(end- t0_p1)/n.int)
      root        <- try(uniroot(fun, interval = interval)$root, silent = TRUE)
      roots_LB[i] <- if(inherits(root, "try-error")) {NA} else root
    }
    roots_LB <-  sort(unique(roots_LB[!is.na(roots_LB)]))
    roots_UB_vec <- sort(unique(c(t0_p1,roots_UB,end)))
    roots_LB_vec <- sort(unique(c(t0_p1,roots_LB,end)))
    roots_vec    <- sort(unique(c(t0_p1,end,roots_UB, roots_LB)))
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

    fun_stat_sig_vec <- function(x) {if(fun_stat_sig_UB_vec(x) != fun_stat_sig_LB_vec(x)) 0 else fun_stat_sig_UB_vec(x)}
    fun_stat_sig_vec <- Vectorize(fun_stat_sig_vec)

    # draw plots
    x_vals  <- seq(start, end, length.out=7*length(timeVec))

    y_range     <- range(c(betahat_splinefun(x_vals), scb_ub_splinefun(x_vals), scb_lb_splinefun(x_vals)), na.rm=TRUE)
    y_range     <- c(y_range[1] - diff(y_range)*0.12, y_range[2] + diff(y_range)*0.08)
    y_range_fun <- function(x) y_range[1]+(y_range[2]-y_range[1])*(x-start)/(end-start)
    curve(y_range_fun, from=start, to=end, lty=1,  xlab = "Event Time", ylab = "", col="white", n=500,...) #ylim =y_range,

    ShadeBetween <- function(x1, x2, f1, f2, ...) {polygon(c(x1, rev(x2)), c(f1, rev(f2)), ...)}

    polygon(
      x = c(seq(t0_p1, end, length.out = 200), rev(seq(t0_p1, end, length.out = 200))),
      y = c(scb_ub_splinefun(seq(t0_p1, end, length.out = 200)), rev(scb_lb_splinefun(seq(t0_p1, end, length.out = 200)))),
      col =  adjustcolor("#48C9B0", alpha.f = 1),   # transparency
      border = NA
    )

    polygon(
      x = c(seq(start, t0, length.out = 200), rev(seq(start, t0, length.out = 200))),
      y = c(scb_ub_splinefun(seq(start, t0, length.out = 200)), rev(scb_lb_splinefun(seq(start, t0, length.out = 200)))),
      col = adjustcolor("#48C9B0", alpha.f = 0.7),   # transparency
      border = NA
    )

    curve(betahat_splinefun, from=start, to=end, lty=1, add=TRUE, xlab = "Event Time", ylab = "", col="blue",n=500,...)
    points(timeVec[which(timeVec<=t0)], betahat[which(timeVec<=t0)], pch = 16, col = rgb(0,0,1,alpha=0.4), cex=1)
    points(timeVec[which(timeVec>=t0_p1)], betahat[which(timeVec>=t0_p1)], pch = 16, col = "blue", cex=1)
    abline(v=(t0+timeVec[which(timeVec==t0)+1])/2, lty=3, col="blue4", lwd=4)


    if(isTRUE(ci.pre)) {arrows(timeVec[which(timeVec<=t0)], ci_lower[which(timeVec<=t0)], timeVec[which(timeVec<=t0)], ci_upper[which(timeVec<=t0)], angle = 90, code = 3, length = 0.025, col = rgb(0,0,1,alpha=0.4));}
    if(isTRUE(ci.post)) {arrows(timeVec[which(timeVec>=t0_p1)], ci_lower[which(timeVec>=t0_p1)], timeVec[which(timeVec>=t0_p1)], ci_upper[which(timeVec>=t0_p1)], angle = 90, code = 3, length = 0.025, col = "blue4");}

    if(!isTRUE(ref.band.pre)) {
      segments(x0=t0_p1, y0=0, x1=end, y1=0, lty=3, lwd=4, col="red")
      ShadeBetween(timeVec[timeVec>=t0_p1], timeVec[timeVec>=t0_p1], honest_ub_splinefun(timeVec[timeVec>=t0_p1]), honest_lb_splinefun(timeVec[timeVec>=t0_p1]), col=rgb(1,0,0,alpha=0.4), border=NA)
      roots_vec <- sort(unique(c(roots_vec, timeVec[which(timeVec==t0_p1)])))
      for (i in 1:(length(roots_vec)-1) ) {
        if ( fun_stat_sig_vec((roots_vec[i]+roots_vec[i+1])/2)==1 & roots_vec[i]>=t0_p1 ) {
          if(isTRUE(verbose)) cat("Uniformly and Positively Significant Time Span: [", roots_vec[i], ",", roots_vec[i+1], "]\n")
        } else {
          if (fun_stat_sig_vec((roots_vec[i]+roots_vec[i+1])/2)==-1 & roots_vec[i]>=t0_p1 ) {
            if(isTRUE(verbose)) cat("Uniformly and Negatively Significant Time Span: [", roots_vec[i], ",", roots_vec[i+1], "]\n")
          }
        }
      }
    } else {
      segments(x0=start, y0=0, x1=end, y1=0, lty=3, lwd=4, col="red")
      ShadeBetween(timeVec[timeVec>=t0_p1], timeVec[timeVec>=t0_p1], honest_ub_splinefun(timeVec[timeVec>=t0_p1]), honest_lb_splinefun(timeVec[timeVec>=t0_p1]), col=rgb(1,0,0,alpha=0.4), border=NA)
      ShadeBetween(timeVec[timeVec<=ifelse(is.null(ta.ts), t0, ta.ts)], timeVec[timeVec<=ifelse(is.null(ta.ts), t0, ta.ts)], honest_ub_splinefun(timeVec[timeVec<=ifelse(is.null(ta.ts), t0, ta.ts)]), honest_lb_splinefun(timeVec[timeVec<=ifelse(is.null(ta.ts), t0, ta.ts)]), col=rgb(1,0,0,alpha=0.2), border=NA)

      roots_vec <- sort(unique(c(roots_vec, timeVec[which(timeVec==t0_p1)])))
      roots_vec <- roots_vec[roots_vec>=timeVec[which(timeVec==t0_p1)]]

      for (i in 1:(length(roots_vec)-1) ) {
        if ( fun_stat_sig_vec((roots_vec[i]+roots_vec[i+1])/2)==1 & roots_vec[i]>=t0_p1) {
          if(isTRUE(verbose)) cat("Uniformly and Positively Significant Time Span: [", roots_vec[i], ",", roots_vec[i+1], "]\n")
        } else {
          if (fun_stat_sig_vec((roots_vec[i]+roots_vec[i+1])/2)==-1 & roots_vec[i]>=t0_p1) {
            if(isTRUE(verbose)) cat("Uniformly and Negatively Significant Time Span: [", roots_vec[i], ",", roots_vec[i+1], "]\n")
          }
        }
      }
    }

    # write a legend for the plot
    if(!is.null(pos.legend)) {

      xlim <- par("usr")[1:2]
      ylim <- par("usr")[3:4]

      upper_text <- "Honest Reference Band"; lower_text <- "Classical Reference Line          "

      text_width_upper   <- strwidth(upper_text, cex=scale.legend) # calculate the width of the text
      text_width_lower   <- strwidth(lower_text, cex=scale.legend)
      fixed_rect_b_width <- strwidth("Classical", cex=scale.legend) # define a longer fixed width for the red rectangle and red dashed line
      fixed_line_width   <- strwidth("Classical", cex=scale.legend)
      extra_padding      <- strwidth("   ", cex=scale.legend) # add extra padding for the white rectangle to avoid text touching borders
      rect_width         <- max(text_width_upper, text_width_lower)+extra_padding+fixed_rect_b_width  # calculate the width of the white rectangle
      rect_height        <- diff(ylim)*0.1*scale.legend # define the height for the white rectangle
      rect_x1            <- mean(xlim)-rect_width/2 # define positions for the white rectangle
      rect_x2            <- mean(xlim)+rect_width/2

      if(pos.legend=="top")    {rect_y2 <- ylim[2]-rect_height*0.5; rect_y1 <- rect_y2-rect_height}
      if(pos.legend=="bottom") {rect_y2 <- ylim[1]+rect_height*0.5; rect_y1 <- rect_y2+rect_height}

      rect(rect_x1, rect_y1, rect_x2, rect_y2, col = "white", border = "black") # draw white rectangular with black border
      left_padding <- rect_width*0.05 # define padding to prevent touching the left border
      rect_b_x1    <- rect_x1+left_padding  # position and draw the fixed-width red rectangle
      rect_b_x2    <- rect_b_x1+fixed_rect_b_width

      if(pos.legend=="top")    {rect_b_y2 <- rect_y2-rect_height*0.15; rect_b_y1 <- rect_b_y2-rect_height*0.35}
      if(pos.legend=="bottom") {rect_b_y2 <- rect_y1-rect_height*0.5; rect_b_y1 <- rect_b_y2+rect_height*0.35}

      if (isTRUE(ref.band.pre)) {
        rect_b_xmid <- (rect_b_x1 + rect_b_x2) / 2   # compute the middle coordinate
        rect(rect_b_x1, rect_b_y1, rect_b_xmid, rect_b_y2, col = rgb(1, 0, 0, alpha = 0.2),  border = NA) # draw the left half (more transparent red)
        rect(rect_b_xmid, rect_b_y1, rect_b_x2, rect_b_y2, col = rgb(1, 0, 0, alpha = 0.4), border = NA)  # draw the right half (less transparent red)
      } else{
        rect(rect_b_x1, rect_b_y1, rect_b_x2, rect_b_y2, col = rgb(1, 0, 0, alpha = 0.4), border = NA)  # draw red rectangle (transparent, no border)
      }

      text(rect_b_x2 + left_padding, (rect_b_y1 + rect_b_y2) / 2, upper_text, cex = scale.legend,  adj = 0)  # add upper text next to red rectangle
      line_x1 <- rect_b_x1  # define position for the red dashed line (aligned with lower text) # same as red rectangle start
      line_x2 <- rect_b_x2  # same fixed width as red rectangle

      if(pos.legend=="top")    {line_y <- rect_y1+rect_height*0.3}
      if(pos.legend=="bottom") {line_y <- rect_y2+rect_height*0.3}

      segments(line_x1, line_y, line_x2, line_y, col = "red", lwd = 4*scale.legend, lty = 3)  # draw the red dashed line with fixed width
      text(line_x2 + left_padding, line_y, lower_text, cex = scale.legend, adj = 0)  # add lower text next to the red dashed line
    }

    # write note on top of the plot
    if (isTRUE(note.pre) & !isTRUE(note.post)) {

      x_left   <- par("usr")[1]
      x_right  <- par("usr")[2]
      y_top    <- par("usr")[4]
      y_bottom <- par("usr")[3]

      box_height <- 0.1 * (y_top - y_bottom)
      x_mid      <- ifelse(is.null(ta.ts), t0, ta.ts)

      par(xpd = TRUE) # Turn off clipping so we can draw *outside* the plot
      segments(x_mid,  y_top,  x_mid,  y_top + box_height, col = "black")

      text((x_left + x_mid)/2, y_top + 0.5*box_height, "Validating", cex = scale.legend)
      par(xpd = FALSE)
    }

    if (!isTRUE(note.pre) & isTRUE(note.post)) {

      x_left   <- par("usr")[1]
      x_right  <- par("usr")[2]
      y_top    <- par("usr")[4]
      y_bottom <- par("usr")[3]

      box_height <- 0.1 * (y_top - y_bottom)
      x_mid      <- timeVec[which(timeVec==t0)+1]

      par(xpd = TRUE)
      segments(x_mid,  y_top,  x_mid,  y_top + box_height, col = "black")
      text((x_mid + x_right)/2, y_top + 0.5*box_height, "Testing", cex = scale.legend)

      par(xpd = FALSE)
    }

    if (isTRUE(note.pre) & isTRUE(note.post)) {

      x_left   <- par("usr")[1]
      x_right  <- par("usr")[2]
      y_top    <- par("usr")[4]
      y_bottom <- par("usr")[3]

      box_height <- 0.1 * (y_top - y_bottom)
      x_mid_pre  <-  ifelse(is.null(ta.ts), t0, ta.ts)
      x_mid_post <- timeVec[which(timeVec==t0)+1]

      par(xpd = TRUE)
      segments(x_mid_pre,  y_top,  x_mid_pre,  y_top + box_height, col = "black")
      segments(x_mid_post,  y_top,  x_mid_post,  y_top + box_height, col = "black")

      text((x_left + x_mid_pre)/2, y_top + 0.5*box_height, "Validating", cex = scale.legend)
      text((x_mid_post + x_right)/2, y_top + 0.5*box_height, "Testing", cex = scale.legend)

      par(xpd = FALSE)
    }
}
  options(warn=1)
}
