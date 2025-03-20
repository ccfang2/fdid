#' This function computes the estimate of roughness parameter function tau(t) using the covariance function of the functional data
#' @description This function computes the estimate of roughness parameter function tau(t) using the covariance function of the functional data.
#' This is adapted from the \link[ffscb]{cov2tau_fun} in package \pkg{ffscb}, but the range of grid is not limited to [0,1].
#'
#' @param cov_mat matrix (pxp) of evaluated covariance function (p=number of discretization point).
#' @param grid.min min value of grid.
#' @param grid.max max value of grid.
#' @param warn option for printing warnings.
#'
#' @return estimate of the roughness parameter function tau(t).
#' @export
#'
#' @references Liebl, D. and M. Reimherr (2023). Fast and fair simultaneous confidence bands for functional parameters. Journal of
#' the Royal Statistical Society Series B: Statistical Methodology 85(3), 842–868
#'
#' @seealso \link[ffscb]{cov2tau_fun}
cov2tau_fun2 <- function(cov_mat, grid.min, grid.max, warn = FALSE){
  ## 'cov_mat' denotes the sample covariance function computed from the sample functions X_1(t),...X_n(t)
  ## Caution: assumed grid is in [0,1]
  ##
  p        <- ncol(cov_mat)
  grid     <- seq(grid.min, grid.max, len=p)
  corr_mat <- stats::cov2cor(cov_mat)
  ## computing the numeric approximation to c_12(t,s) := \partial^2 c(t,s)/(\partial t \partial t)
  ## with c_12(t,s) evaluated at t=t and s=t
  a1  <- 1 #Equals 1: corr_mat[cbind(2:p,2:p)]        # corr(t+h, t+h) with h=0.005, and t\in{0.005,0.015,0.025,...,0.995}
  a2  <- corr_mat[cbind(1:(p-1),2:p)]                 # corr(t-h, t+h)
  a3  <- 1 #Equals 1: corr_mat[cbind(1:(p-1),1:(p-1))]# corr(t-h, t-h)
  ##
  h   <- diff(grid)[1]/2
  suppressWarnings(
    tau <- sqrt( c(a1 - 2*a2 + a3) /( 4* h^2 ) )
  )
  if(any(is.na(tau))){
    if(warn){
      warning("\nCaution there were negative tau-values!\n
             Probably, this was caused by a cov-estimate\n
             computed from fragmentary functions.\n
             We use only the positive values and fill the \n
             missing tau-values by approximations using linear interpolations.")
    }
    if(is.na(tau[1]  )){tau[1]   <- utils::head(c(stats::na.omit(tau)),n = 1)}
    if(is.na(tau[p-1])){tau[p-1] <- utils::tail(c(stats::na.omit(tau)),n = 1)}
    tau <- stats::approx(x=c(1:(p-1))[!is.na(tau)], y=tau[!is.na(tau)], xout = c(1:(p-1)))$y
  }
  ## tau has length p-1, so we interpolate to get length p
  xx    <- c(grid.min,      grid[-p] + h, grid.max               )
  yy    <- c(tau[1], tau,          tau[length(tau)])
  tau_t <- stats::approx(y=yy, x=xx, xout = grid)$y
  ##
  return(tau_t)
}


