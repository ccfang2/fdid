#' Custom Plot Method for objects of class \code{"fdid_scb"}
#' @description This function is an S3 method for \link[base]{plot}, specifically designed for objects of class \code{"fdid_scb"}.
#' It plots simultaneous confidence band with bounds for honest inference.
#'
#' @param object an object of class \code{"fdid_scb"}. The object to be plotted.
#' @param ta.ts a numeric value that indicates the time point based on which the sample estimates for treatment anticipation biases are computed. It should be less or equal to the time point after which units start responding to a future treatment. It is supposed to be greater than the minimal event time and no greater than t0 in \code{object}. If it is NULL, there is no anticipation.
#' @param ta.s a numeric value of control parameter for honest inference under treatment anticipation. If NULL, the critical value \eqn{t_{1-\alpha/2, df}} or \eqn{z_{1-\alpha/2}} is applied with the same \eqn{\alpha} used for constructing the pointwise confidence intervals. If no treatment anticipation is defined in \code{ta.ts}, the value of \code{ta.s} is ignored.
#' @param frm.mbar a numeric value of control parameter for honest inference under functional relative magnitudes.
#' @param ftr.m a numeric value of control parameter for honest inference under functional trend restrictions.
#' @param frmtr.mbar a numeric value of control parameter for combining functional relative magnitudes and functional trend restrictions.
#' @param ref.band.pre a logical value. If TRUE, the reference band for pre-treatment periods is also plotted.
#' @param note.pre a logical value. If TRUE, the note for pre-treatment periods is given on top of plot.
#' @param note.post a logical value. If TRUE, the note for post-treatment periods is given on top of plot.
#' @param ci a logical value. If TRUE, the point-wise confidence intervals are also plotted.
#' @param pos.legend a character value of "top" or "bottom" that indicates the position of legend. If NULL, the legend is not printed.
#' @param scale.legend a positive number that defines the size of legend. If \code{pos.legend} is NULL, the value of \code{scale.legend} is ignored.
#' @param ... Additional arguments to be passed to \link[base]{plot}.
#'
#' @return The function returns a plot with simultaneous confidence band for event study coefficients in a functional framework, together with bounds for honest inference, if properly defined.
#'
#' @import pracma
#' @export
#'
#' @references Fang, C. and Liebl, D. (2025). Making Event Study Plots Honest: A Functional Data Approach to Causal Inference.
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
#' ## adding a label for Y-axis
#' title(ylab="Effects of Duty-to-Bargain Laws")
#'
#' ## honest inference under treatment anticipation
#' plot(fdid_scb_est, ta.ts=-3, scale.legend=1.4)
#'
#' ## honest inference under violation of parallel trends
#' plot(fdid_scb_est, frm.mbar=0.4, scale.legend=1.4)
plot.fdid_scb <- function(object,
                          ta.ts=NULL,
                          ta.s=NULL,
                          frm.mbar=NULL,
                          ftr.m=NULL,
                          frmtr.mbar=NULL,
                          ref.band.pre=TRUE,
                          note.pre=TRUE,
                          note.post=TRUE,
                          ci=TRUE,
                          pos.legend="top",
                          scale.legend=1,
                          ...) {

  # give out a warn
  # warning_note <- paste("The inference result around the reference time t=", object$data$t0, "(i.e. between two event time closest to the reference time) should be treated with caution.")
  # warning(warning_note)

  options(warn=-1)

  # check conditions
  if (!base::inherits(object,"fdid_scb")) stop("The input 'object' should be an output of function 'fdid_scb'.")
  if (!is.null(ta.ts) && (!is.numeric(ta.ts) || length(ta.ts) != 1)) stop("The input 'ta.ts' should be either NULL or a numeric scalar.")
  if (!is.null(ta.ts) && !(ta.ts %in% object$scb$event_t[which(object$scb$event_t<=object$data$t0)])) stop("If not NULL, the input 'ta.ts' should be among the pre-treatment event time in 'object'.")
  if (!is.null(ta.ts) && ta.ts == object$scb$event_t[1] && !all(sapply(list(frm.mbar, ftr.m, frmtr.mbar), is.null))) stop("If 'ta.ts' is defined to be the first event time, 'frm.mbar', 'ftr.m' and 'frmtr.mbar' must be NULL, because there is no available data for computing pre-trend differences.")
  if (!is.null(ta.s) && (!is.numeric(ta.s) || length(ta.s) != 1 || ta.s < 0)) stop("The input 'ta.s' should be either NULL or a numeric non-negative scalar.")
  if (!is.null(frm.mbar) && (!is.numeric(frm.mbar) || length(frm.mbar) != 1 || frm.mbar < 0)) stop("The input 'frm.mbar' should be either NULL or a numeric non-negative scalar.")
  if (!is.null(ftr.m) && (!is.numeric(ftr.m) || length(ftr.m) != 1 || ftr.m < 0)) stop("The input 'ftr.m' should be either NULL or a numeric non-negative scalar.")
  if (!is.null(frmtr.mbar) && (!is.numeric(frmtr.mbar) || length(frmtr.mbar) != 1 || frmtr.mbar < 0)) stop("The input 'frmtr.mbar' should be either NULL or a numeric non-negative scalar.")
  if (!is.null(frm.mbar) && !is.null(ftr.m) && !is.null(frmtr.mbar)) stop("Exactly one of inputs 'frm.mbar', 'ftr.m' and 'frmtr.mbar' can have a value, or all three are NULL.")
  if (!is.logical(ref.band.pre)) stop("The input 'ref.band.pre' should be logical.")
  if (!is.logical(note.pre)) stop("The input 'note.pre' should be logical.")
  if (!is.logical(note.post)) stop("The input 'note.post' should be logical.")
  if (!is.logical(ci)) stop("The input 'ci' should be logical.")
  if (!is.null(pos.legend) && !pos.legend %in% c("top", "bottom")) stop("The input 'pos.legend' must be 'top', 'bottom', or NULL.")
  if (!is.numeric(scale.legend) || length(scale.legend) != 1 || scale.legend <= 0 || !is.finite(scale.legend)) stop("The input 'scale.legend' must be a positive number.")

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

  if (is.null(ta.ts)) ta.ts <- t0

  # regular statistical inference
  if (ta.ts==t0 & is.null(frm.mbar) & is.null(ftr.m) & is.null(frmtr.mbar)) {

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
    y_range <- c(y_range[1] - diff(y_range)*0.08, y_range[2] + diff(y_range)*0.08)

    # draw plot
    plot(timeVec, betahat,  type = "p", ylim =y_range, xlab = "Event Time", ylab = "", pch = 16, col = "blue", ...)
    if(isTRUE(ci)) arrows(timeVec, ci_lower, timeVec, ci_upper, angle = 90, code = 3, length = 0.025, col = "blue4")
    #grid()

    curve(betahat_splinefun, add=TRUE, lty=1)
    curve(scb_ub_splinefun, add=TRUE, lty=1)
    curve(scb_lb_splinefun, add=TRUE, lty=1)

    if(!isTRUE(ref.band.pre)) {
      segments(x0=t0, y0=0, x1=end, y1=0, lty=3, lwd=4, col="red")
      roots_vec_temp <- sort(unique(c(roots_vec, timeVec[which(timeVec==t0)-1], timeVec[which(timeVec==t0)+1])))

      intervals <- data.frame(
        left  = head(roots_vec_temp, -1),
        right = tail(roots_vec_temp, -1)
      )

      stat_sig_vec_temp <- rep(NA, nrow(intervals))

      for (i in 1:nrow(intervals)) {
        l <- intervals$left[i]
        r <- intervals$right[i]
        if (l >= timeVec[which(timeVec==t0)-1] && r <= timeVec[which(timeVec==t0)+1]) {
          stat_sig_vec_temp[i] <- 0
        } else {
          idx <- which(roots_vec <= l)[length(which(roots_vec  <= l))]
          stat_sig_vec_temp[i] <- stat_sig_vec[idx]
        }
      }

      stat_sig_vec <- stat_sig_vec_temp
      roots_vec    <- roots_vec_temp

      # if(!(t0 %in% roots_vec)) {
      #   roots_vec <- sort(c(roots_vec, t0))
      #   stat_sig_vec <- append(stat_sig_vec, stat_sig_vec[which(roots_vec==t0)-1], after=which(roots_vec==t0)-1-1)
      # }

      for (i in 1:(length(roots_vec)-1)) {
        if (stat_sig_vec[i]==1 & roots_vec[i]>=t0 ) {
          rect(roots_vec[i], y_range[1], roots_vec[i+1], y_range[2], col=rgb(0.5, 0.5, 0.5, alpha=0.4), border=NA) #alpha=0.7
        } else {
          if (stat_sig_vec[i]==-1 & roots_vec[i]>=t0) {
            rect(roots_vec[i], y_range[1], roots_vec[i+1], y_range[2], col=rgb(0.5, 0.5, 0.5, alpha=0.4), border=NA)
          }
        }
      }
    } else {
      segments(x0=start, y0=0, x1=end, y1=0, lty=3, lwd=4, col="red")

      roots_vec_temp <- sort(unique(c(roots_vec, timeVec[which(timeVec==t0)-1], timeVec[which(timeVec==t0)+1])))

      intervals <- data.frame(
        left  = head(roots_vec_temp, -1),
        right = tail(roots_vec_temp, -1)
      )

      stat_sig_vec_temp <- rep(NA, nrow(intervals))

      for (i in 1:nrow(intervals)) {
        l <- intervals$left[i]
        r <- intervals$right[i]
        if (l >= timeVec[which(timeVec==t0)-1] && r <= timeVec[which(timeVec==t0)+1]) {
          stat_sig_vec_temp[i] <- 0
        } else {
          idx <- which(roots_vec <= l)[length(which(roots_vec  <= l))]
          stat_sig_vec_temp[i] <- stat_sig_vec[idx]
        }
      }

      stat_sig_vec <- stat_sig_vec_temp
      roots_vec    <- roots_vec_temp

      for (i in 1:(length(roots_vec)-1)) {
        if (stat_sig_vec[i]==1) {
          rect(roots_vec[i], y_range[1], roots_vec[i+1], y_range[2], col=rgb(0.5, 0.5, 0.5, alpha=0.4), border=NA) #alpha=0.7
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

      #text              <- "Reference value          "
      text              <- "Classical Reference Line          "
      text_width        <- strwidth(text, cex = scale.legend)
      #fixed_line_width  <- strwidth("Reference", cex=scale.legend)
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

    # write note on top of the plot
    if (isTRUE(note.pre) & !isTRUE(note.post)) {

      # Get the user coordinates
      x_left <- par("usr")[1]
      x_right <- par("usr")[2]
      y_top <- par("usr")[4]
      y_bottom <- par("usr")[3]

      box_height <- 0.1 * (y_top - y_bottom)  # height = 10% of plot range
      x_mid      <- t0                        # border between boxes

      par(xpd = TRUE) # Turn off clipping so we can draw *outside* the plot
      segments(x_mid,  y_top,  x_mid,  y_top + box_height, col = "black")
      text((x_left + x_mid)/2, y_top + 0.5*box_height, "Validating", cex = scale.legend)
      par(xpd = FALSE) # Turn clipping back on
    }

    if (!isTRUE(note.pre) & isTRUE(note.post)) {

      # Get the user coordinates
      x_left <- par("usr")[1]
      x_right <- par("usr")[2]
      y_top <- par("usr")[4]
      y_bottom <- par("usr")[3]

      box_height <- 0.1 * (y_top - y_bottom)  # height = 10% of plot range
      x_mid      <- t0                         # border between boxes

      par(xpd = TRUE) # Turn off clipping so we can draw *outside* the plot
      segments(x_mid,  y_top,  x_mid,  y_top + box_height, col = "black")
      text((x_mid + x_right)/2, y_top + 0.5*box_height, "Testing", cex = scale.legend)
      par(xpd = FALSE) # Turn clipping back on
    }

    if (isTRUE(note.pre) & isTRUE(note.post)) {

      # Get the user coordinates
      x_left <- par("usr")[1]
      x_right <- par("usr")[2]
      y_top <- par("usr")[4]
      y_bottom <- par("usr")[3]

      box_height <- 0.1 * (y_top - y_bottom)  # height = 10% of plot range
      x_mid      <- t0                         # border between boxes

      par(xpd = TRUE) # Turn off clipping so we can draw *outside* the plot
      segments(x_mid,  y_top,  x_mid,  y_top + box_height, col = "black")
      text((x_left + x_mid)/2, y_top + 0.5*box_height, "Validating", cex = scale.legend)
      text((x_mid + x_right)/2, y_top + 0.5*box_height, "Testing", cex = scale.legend)
      par(xpd = FALSE) # Turn clipping back on
    }

  } else {

    # honest inference
    if(is.null(frm.mbar) & is.null(ftr.m) & is.null(frmtr.mbar)) {

      if(is.null(ta.s)) {
        ta_ub <- ci_upper[names(ci_upper)==as.character(ta.ts)]
        ta_lb <- ci_lower[names(ci_lower)==as.character(ta.ts)]
      } else {
        ta_ub <- betahat[which(timeVec==ta.ts)]+ta.s*sqrt(diag(covhat)[which(timeVec==ta.ts)])
        ta_lb <- betahat[which(timeVec==ta.ts)]-ta.s*sqrt(diag(covhat)[which(timeVec==ta.ts)])
      }

      honest_ub_splinefun <- function(x) ta_ub
      honest_lb_splinefun <- function(x) ta_lb

      honest_ub_splinefun <- Vectorize(honest_ub_splinefun)
      honest_lb_splinefun <- Vectorize(honest_lb_splinefun)

    }


    if(!is.null(frm.mbar) & is.null(ftr.m) & is.null(frmtr.mbar)) {

      if(is.null(ta.s)) {
        ta_ub <- ci_upper[names(ci_upper)==as.character(ta.ts)]
        ta_lb <- ci_lower[names(ci_lower)==as.character(ta.ts)]
      } else {
        ta_ub <- betahat[which(timeVec==ta.ts)]+ta.s*sqrt(diag(covhat)[which(timeVec==ta.ts)])
        ta_lb <- betahat[which(timeVec==ta.ts)]-ta.s*sqrt(diag(covhat)[which(timeVec==ta.ts)])
      }

      # frm bounds
      ts             <- seq(start, ta.ts, length.out=50*sum(timeVec<=ta.ts))
      betas_deriv    <- pracma::fderiv(betahat_splinefun, ts, n=1, method="backward")
      mean_abs_deriv <- mean(abs(betas_deriv))

      frm_ub_ta_ub <- function (x) {ifelse(x>=ta.ts, frm.mbar*mean_abs_deriv*(x-ta.ts)+ta_ub, -frm.mbar*mean_abs_deriv*(x-ta.ts)+ta_ub)}
      frm_lb_ta_ub <- function (x) {ifelse(x>=ta.ts, -frm.mbar*mean_abs_deriv*(x-ta.ts)+ta_ub, frm.mbar*mean_abs_deriv*(x-ta.ts)+ta_ub)}
      frm_ub_ta_lb <- function (x) {ifelse(x>=ta.ts, frm.mbar*mean_abs_deriv*(x-ta.ts)+ta_lb, -frm.mbar*mean_abs_deriv*(x-ta.ts)+ta_lb)}
      frm_lb_ta_lb <- function (x) {ifelse(x>=ta.ts, -frm.mbar*mean_abs_deriv*(x-ta.ts)+ta_lb, frm.mbar*mean_abs_deriv*(x-ta.ts)+ta_lb)}

      x_vals <- seq(start, end, length.out = 5*length(timeVec))

      if (ta.ts!=t0) {
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
        ta_ub <- ci_upper[names(ci_upper)==as.character(ta.ts)]
        ta_lb <- ci_lower[names(ci_lower)==as.character(ta.ts)]
      } else {
        ta_ub <- betahat[which(timeVec==ta.ts)]+ta.s*sqrt(diag(covhat)[which(timeVec==ta.ts)])
        ta_lb <- betahat[which(timeVec==ta.ts)]-ta.s*sqrt(diag(covhat)[which(timeVec==ta.ts)])
      }

      # ftr bounds
      ts          <- seq(start, ta.ts, length.out=50*sum(timeVec<=ta.ts))
      betas_deriv <- pracma::fderiv(betahat_splinefun, ts, n=1, method="backward")
      slope       <- mean(betas_deriv)

      ftr_ub_ta_ub <- function (x) {ifelse(x>=ta.ts, (ftr.m*abs(slope)+slope)*(x-ta.ts)+ta_ub, (-ftr.m*abs(slope)+slope)*(x-ta.ts)+ta_ub)}
      ftr_lb_ta_ub <- function (x) {ifelse(x>=ta.ts, (-ftr.m*abs(slope)+slope)*(x-ta.ts)+ta_ub, (ftr.m*abs(slope)+slope)*(x-ta.ts)+ta_ub)}
      ftr_ub_ta_lb <- function (x) {ifelse(x>=ta.ts, (ftr.m*abs(slope)+slope)*(x-ta.ts)+ta_lb, (-ftr.m*abs(slope)+slope)*(x-ta.ts)+ta_lb)}
      ftr_lb_ta_lb <- function (x) {ifelse(x>=ta.ts, (-ftr.m*abs(slope)+slope)*(x-ta.ts)+ta_lb, (ftr.m*abs(slope)+slope)*(x-ta.ts)+ta_lb)}

      x_vals <- seq(start, end, length.out=5*length(timeVec))

      if (ta.ts!=t0) {
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
        ta_ub <- ci_upper[names(ci_upper)==as.character(ta.ts)]
        ta_lb <- ci_lower[names(ci_lower)==as.character(ta.ts)]
      } else {
        ta_ub <- betahat[which(timeVec==ta.ts)]+ta.s*sqrt(diag(covhat)[which(timeVec==ta.ts)])
        ta_lb <- betahat[which(timeVec==ta.ts)]-ta.s*sqrt(diag(covhat)[which(timeVec==ta.ts)])
      }

      # frmtr bounds
      ts             <- seq(start, ta.ts, length.out=50*sum(timeVec<=ta.ts))
      betas_deriv    <- pracma::fderiv(betahat_splinefun, ts, n=1, method="backward")
      mean_abs_deriv <- mean(abs(betas_deriv))
      slope          <- mean(betas_deriv)

      ftr_ub_ta_ub <- function (x) {ifelse(x>=ta.ts, (frmtr.mbar*mean_abs_deriv+slope)*(x-ta.ts)+ta_ub, (-frmtr.mbar*mean_abs_deriv+slope)*(x-ta.ts)+ta_ub)}
      ftr_lb_ta_ub <- function (x) {ifelse(x>=ta.ts, (-frmtr.mbar*mean_abs_deriv+slope)*(x-ta.ts)+ta_ub, (frmtr.mbar*mean_abs_deriv+slope)*(x-ta.ts)+ta_ub)}
      ftr_ub_ta_lb <- function (x) {ifelse(x>=ta.ts, (frmtr.mbar*mean_abs_deriv+slope)*(x-ta.ts)+ta_lb, (-frmtr.mbar*mean_abs_deriv+slope)*(x-ta.ts)+ta_lb)}
      ftr_lb_ta_lb <- function (x) {ifelse(x>=ta.ts, (-frmtr.mbar*mean_abs_deriv+slope)*(x-ta.ts)+ta_lb, (frmtr.mbar*mean_abs_deriv+slope)*(x-ta.ts)+ta_lb)}

      x_vals <- seq(start, end, length.out=5*length(timeVec))

      if (ta.ts!=t0) {
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

    fun_stat_sig_vec <- function(x) {if(fun_stat_sig_UB_vec(x) != fun_stat_sig_LB_vec(x) || (x >= timeVec[which(timeVec==t0)-1] && x <= timeVec[which(timeVec==t0)+1])) 0 else fun_stat_sig_UB_vec(x)}
    fun_stat_sig_vec <- Vectorize(fun_stat_sig_vec)

    # draw plots
    x_vals  <- seq(start, end, length.out=5*length(timeVec))
    y_range <- range(c(betahat_splinefun(x_vals), scb_ub_splinefun(x_vals), scb_lb_splinefun(x_vals)), na.rm=TRUE)
    y_range <- c(y_range[1] - diff(y_range)*0.08, y_range[2] + diff(y_range)*0.08)

    plot(timeVec, betahat,  type = "p", ylim =y_range, xlab = "Event Time", ylab = "", pch = 16, col = "blue", ...)
    #grid()

    if(isTRUE(ci)) arrows(timeVec, ci_lower, timeVec, ci_upper, angle = 90, code = 3, length = 0.025, col = "blue4")
    curve(betahat_splinefun, add=TRUE, lty=1)
    curve(scb_ub_splinefun, add=TRUE, lty=1)
    curve(scb_lb_splinefun, add=TRUE, lty=1)

    ShadeBetween <- function(x1, x2, f1, f2, ...) {polygon(c(x1, rev(x2)), c(f1, rev(f2)), ...)}

    if(!isTRUE(ref.band.pre)) {
      segments(x0=t0, y0=0, x1=end, y1=0, lty=3, lwd=4, col="red")
      ShadeBetween(timeVec[timeVec>=t0], timeVec[timeVec>=t0], honest_ub_splinefun(timeVec[timeVec>=t0]), honest_lb_splinefun(timeVec[timeVec>=t0]), col=rgb(1,0,0,alpha=0.4), border=NA)

      roots_vec <- sort(unique(c(roots_vec, timeVec[which(timeVec==t0)-1], timeVec[which(timeVec==t0)+1])))
      # if(!(t0 %in%  roots_vec)) {roots_vec <- sort(c(roots_vec, t0))}

      for (i in 1:(length(roots_vec)-1) ) {
        if ( fun_stat_sig_vec((roots_vec[i]+roots_vec[i+1])/2)==1 & roots_vec[i]>=t0 ) {
          rect(roots_vec[i], y_range[1], roots_vec[i+1], y_range[2], col=rgb(0.5, 0.5, 0.5, alpha=0.4), border=NA) #alpha=0.7
        } else {
          if (fun_stat_sig_vec((roots_vec[i]+roots_vec[i+1])/2)==-1 & roots_vec[i]>=t0 ) {
            rect(roots_vec[i], y_range[1], roots_vec[i+1], y_range[2], col=rgb(0.5, 0.5, 0.5, alpha=0.4), border=NA)
          }
        }
      }
    } else {

      segments(x0=start, y0=0, x1=end, y1=0, lty=3, lwd=4, col="red")
      ShadeBetween(timeVec[timeVec>=t0], timeVec[timeVec>=t0], honest_ub_splinefun(timeVec[timeVec>=t0]), honest_lb_splinefun(timeVec[timeVec>=t0]), col=rgb(1,0,0,alpha=0.4), border=NA)
      ShadeBetween(timeVec[timeVec<=t0], timeVec[timeVec<=t0], honest_ub_splinefun(timeVec[timeVec<=t0]), honest_lb_splinefun(timeVec[timeVec<=t0]), col=rgb(1,0,0,alpha=0.2), border=NA)

      roots_vec <- sort(unique(c(roots_vec, timeVec[which(timeVec==t0)-1], timeVec[which(timeVec==t0)+1])))
      roots_vec <- roots_vec[roots_vec >= timeVec[which(timeVec==t0)+1] ]

      for (i in 1:(length(roots_vec)-1) ) {
        if ( fun_stat_sig_vec((roots_vec[i]+roots_vec[i+1])/2)==1  ) {
          rect(roots_vec[i], y_range[1], roots_vec[i+1], y_range[2], col=rgb(0.5, 0.5, 0.5, alpha=0.4), border=NA) #alpha=0.7
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

      if (ta.ts!=t0 & (!is.null(frm.mbar) | !is.null(ftr.m) | !is.null(frmtr.mbar))) {upper_text <- "Honest Reference Band"; lower_text <- "Classical Reference Line          "}
      if (ta.ts==t0 & (!is.null(frm.mbar) | !is.null(ftr.m) | !is.null(frmtr.mbar))) {upper_text <- "Honest Reference Band"; lower_text <- "Classical Reference Line          "}
      if (ta.ts!=t0 & is.null(frm.mbar) & is.null(ftr.m) & is.null(frmtr.mbar)) {upper_text <- "Honest Reference Band"; lower_text <- "Classical Reference Line          "}

      text_width_upper   <- strwidth(upper_text, cex=scale.legend) # calculate the width of the text dynamically
      text_width_lower   <- strwidth(lower_text, cex=scale.legend)
      #fixed_rect_b_width <- strwidth("Reference", cex=scale.legend) # define a longer fixed width for the red rectangle (B) and red dashed line
      fixed_rect_b_width <- strwidth("Classical", cex=scale.legend) # define a longer fixed width for the red rectangle (B) and red dashed line
      #fixed_line_width   <- strwidth("Reference", cex=scale.legend)
      fixed_line_width   <- strwidth("Classical", cex=scale.legend)
      extra_padding      <- strwidth("   ", cex=scale.legend) # add extra padding for the white rectangle to avoid text touching borders
      rect_width         <- max(text_width_upper, text_width_lower)+extra_padding+fixed_rect_b_width  # calculate the width of the white rectangle (A) based on the longest text
      rect_height        <- diff(ylim)*0.1*scale.legend # define the height for the white rectangle (A)
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


      if (isTRUE(ref.band.pre)) {
        rect_b_xmid <- (rect_b_x1 + rect_b_x2) / 2   # compute the middle x coordinate
        rect(rect_b_x1, rect_b_y1, rect_b_xmid, rect_b_y2, col = rgb(1, 0, 0, alpha = 0.2),  border = NA) # draw the left half (more transparent red)
        rect(rect_b_xmid, rect_b_y1, rect_b_x2, rect_b_y2, col = rgb(1, 0, 0, alpha = 0.4), border = NA)  # draw the right half (less transparent red)
      } else{
        rect(rect_b_x1, rect_b_y1, rect_b_x2, rect_b_y2, col = rgb(1, 0, 0, alpha = 0.4), border = NA)  # draw red rectangle B (transparent, no border)
      }

      text(rect_b_x2 + left_padding, (rect_b_y1 + rect_b_y2) / 2, upper_text, cex = scale.legend,  adj = 0)  # add upper text next to red rectangle B

      line_x1 <- rect_b_x1  # define position for the red dashed line (aligned with lower text) # same as red rectangle start
      line_x2 <- rect_b_x2  # same fixed width as red rectangle

      if(pos.legend=="top")    {line_y <- rect_y1+rect_height*0.3}
      if(pos.legend=="bottom") {line_y <- rect_y2+rect_height*0.3}

      segments(line_x1, line_y, line_x2, line_y, col = "red", lwd = 4*scale.legend, lty = 3)  # draw the red dashed line (fixed width)
      text(line_x2 + left_padding, line_y, lower_text, cex = scale.legend, adj = 0)  # add lower text next to the red dashed line
    }

    # write note on top of the plot
    if (isTRUE(note.pre) & !isTRUE(note.post)) {

      # Get the user coordinates
      x_left <- par("usr")[1]
      x_right <- par("usr")[2]
      y_top <- par("usr")[4]
      y_bottom <- par("usr")[3]

      box_height <- 0.1 * (y_top - y_bottom)  # height = 10% of plot range
      x_mid      <- t0                         # border between boxes

      par(xpd = TRUE) # Turn off clipping so we can draw *outside* the plot
      segments(x_mid,  y_top,  x_mid,  y_top + box_height, col = "black")
      text((x_left + x_mid)/2, y_top + 0.5*box_height, "Validating", cex = scale.legend)
      par(xpd = FALSE) # Turn clipping back on
    }

    if (!isTRUE(note.pre) & isTRUE(note.post)) {

      # Get the user coordinates
      x_left <- par("usr")[1]
      x_right <- par("usr")[2]
      y_top <- par("usr")[4]
      y_bottom <- par("usr")[3]

      box_height <- 0.1 * (y_top - y_bottom)  # height = 10% of plot range
      x_mid      <- t0                         # border between boxes

      par(xpd = TRUE) # Turn off clipping so we can draw *outside* the plot
      segments(x_mid,  y_top,  x_mid,  y_top + box_height, col = "black")
      text((x_mid + x_right)/2, y_top + 0.5*box_height, "Testing", cex = scale.legend)
      par(xpd = FALSE) # Turn clipping back on
    }

    if (isTRUE(note.pre) & isTRUE(note.post)) {

      # Get the user coordinates
      x_left <- par("usr")[1]
      x_right <- par("usr")[2]
      y_top <- par("usr")[4]
      y_bottom <- par("usr")[3]

      box_height <- 0.1 * (y_top - y_bottom)  # height = 10% of plot range
      x_mid      <- t0                         # border between boxes

      par(xpd = TRUE) # Turn off clipping so we can draw *outside* the plot
      segments(x_mid,  y_top,  x_mid,  y_top + box_height, col = "black")
      text((x_left + x_mid)/2, y_top + 0.5*box_height, "Validating", cex = scale.legend)
      text((x_mid + x_right)/2, y_top + 0.5*box_height, "Testing", cex = scale.legend)
      par(xpd = FALSE) # Turn clipping back on
    }
  }
  options(warn=1)
}


# Possible Extension: Add an output of significant ranges. This could be done with a generic summary() function. It can
# give out the reference band, the significant regions, and the sign of significance.

