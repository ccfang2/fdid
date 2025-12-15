#' Drawing classical event study plot
#' @description This function is to draw the classical event study plot with point-wise confidence intervals for event study coefficients.
#'
#' @param object an object of class \code{"fdid_scb"}. The object to be plotted.
#' @param pos.legend a character value of "top" or "bottom" that indicates the position of legend. If NULL, the legend is not printed.
#' @param scale.legend a positive number that defines the size of legend. If \code{pos.legend} is NULL, the value of \code{scale.legend} is ignored.
#' @param ... Additional arguments to be passed to \code{EventStudyPlot_Classical}.
#'
#' @return The \code{EventStudyPlot_Classical} function returns a classical event study plot with point-wise confidence intervals for event study coefficients.
#'
#' @export
#'
#' @references Fang, C. and Liebl, D. (2025). Making Event Study Plots Honest: A Functional Data Approach to Causal Inference. \href{https://arxiv.org/abs/2512.06804}{arXiv:2512.06804}.
#' @seealso \link{fdid_scb}
#'
#' @examples
#' data(LWdata)
#' fdid_scb_est <- fdid_scb(beta=LWdata$beta, cov=LWdata$cov, t0=LWdata$t0)
#'
#' par(cex.axis = 1.4, cex.lab = 1.4, cex.main = 1.4)
#' EventStudyPlot_Classical(fdid_scb_est, scale.legend=1.4)
#' title(ylab="Effects of Duty-to-Bargain Laws")
EventStudyPlot_Classical <- function(object,
                                     pos.legend="top",
                                     scale.legend=1,
                                     ...) {

  options(warn=-1)

  # check conditions
  if (!base::inherits(object,"fdid_scb")) stop("The input 'object' should be an output of function 'fdid_scb'.")
  if (!is.null(pos.legend) && !pos.legend %in% c("top", "bottom")) stop("The input 'pos.legend' must be 'top', 'bottom', or NULL.")
  if (!is.numeric(scale.legend) || length(scale.legend) != 1 || scale.legend <= 0 || !is.finite(scale.legend)) stop("The input 'scale.legend' must be a positive number.")

  # extract data from object
  betahat <- object$data$beta[,1]
  covhat  <- object$data$cov
  timeVec <- object$data$beta[,2]

  ci_upper <- object$ci[,"ci_upper"]
  ci_lower <- object$ci[,"ci_lower"]

  start                     <- min(object$data$beta[,2])
  end                       <- max(object$data$beta[,2])
  t0                        <- object$data$t0
  ci.alpha                  <- object$data$ci.alpha

    # draw plot
    y_range     <- range(c(ci_upper, ci_lower), na.rm=TRUE)
    y_range     <- c(y_range[1] - diff(y_range)*0.12, y_range[2] + diff(y_range)*0.08)
    y_range_fun <- function(x) y_range[1]+(y_range[2]-y_range[1])*(x-start)/(end-start)
    curve(y_range_fun, from=start, to=end, lty=1,  xlab = "Event Time", ylab = "", col="white",...)

    points(timeVec, betahat, pch = 16, xlab = "Event Time", ylab = "", col = "blue",...)
    arrows(timeVec, ci_lower, timeVec, ci_upper, angle = 90, code = 3, length = 0.025, col = "blue4")
    segments(x0=timeVec[1],y0=0,x1=timeVec[length(timeVec)],y1=0, lty=3, lwd=4, col="red" )

    # write a legend for the plot
    if(!is.null(pos.legend)) {
        xlim <- par("usr")[1:2]
        ylim <- par("usr")[3:4]

        upper_text <- "Event Study (DiD) Estimates"; lower_text <- "Classical Reference Line            "

        text_width_upper   <- strwidth(upper_text, cex=scale.legend) # calculate the width of the text dynamically
        text_width_lower   <- strwidth(lower_text, cex=scale.legend)
        fixed_rect_b_width <- strwidth("Classical", cex=scale.legend) # define a longer fixed width for the red rectangle (B) and red dashed line
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

        rect_b_xmid <- (rect_b_x1 + rect_b_x2) / 2
        rect_b_ymid <- (rect_b_y1 + rect_b_y2) / 2

        points(rect_b_xmid, rect_b_ymid, pch = 21, bg = "blue", cex = scale.legend)
        text(rect_b_x2 + left_padding, (rect_b_y1 + rect_b_y2) / 2, upper_text, cex = scale.legend,  adj = 0)  # add upper text next to red rectangle B

        line_x1 <- rect_b_x1  # define position for the red dashed line (aligned with lower text) # same as red rectangle start
        line_x2 <- rect_b_x2  # same fixed width as red rectangle

        if(pos.legend=="top")    {line_y <- rect_y1+rect_height*0.3}
        if(pos.legend=="bottom") {line_y <- rect_y2+rect_height*0.3}

        segments(line_x1, line_y, line_x2, line_y, col = "red", lwd = 4*scale.legend, lty = 3)  # draw the red dashed line (fixed width)
        text(line_x2 + left_padding, line_y, lower_text, cex = scale.legend, adj = 0)  # add lower text next to the red dashed line
     }
  options(warn=1)
}
