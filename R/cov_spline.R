#' Natural Spline Interpolation on Covariance
#' @description The \code{cov_spline} function is used to perform natural spline interpolation on the covariance matrix.
#' See Fang and Liebl (2025) for detailed information.
#'
#' @param cov a numeric square matrix of raw covariance. It must not contain missing values.
#' @param grid a numeric vector of grid, at which the raw covariance is calculated.
#' @param n_intrpl the number of points at which the covariance is interpolated. The interpolation is performed on a sequence
#' of \code{n_intrpl} equidistant points from the minimum to the maximum of \code{grid}.
#'
#' @return The \code{cov_spline} function returns a numeric matrix of covariance interpolation.
#' @import splines2
#' @export
#'
#' @references Fang, C. and Liebl, D. (2025). Making Event Study Plots Honest: A Functional Data Approach to Causal Inference.
#'
#' @examples
#' cov_mat <- cov(matrix(rnorm(250),nrow=50))
#' grid <- seq(0,1,len=5)
#' cov_spline_mat <- cov_spline(cov=cov_mat, grid=grid, n_intrpl=50)
cov_spline <- function(cov, grid, n_intrpl) {

  options(warn=-1)

  #check conditions
  if(!is.numeric(cov)) stop("The input 'cov' needs to be numeric.")
  if(ncol(cov)!=nrow(cov)) stop("The input 'cov' needs to be a square matrix.")
  if(!is.numeric(grid)) stop("The input 'grid' needs to be numeric.")
  if(length(grid)!=ncol(cov)) stop("The length of 'grid' needs to be the same as the number of columns/rows of 'cov'.")

  #extract info from grid
  a <- min(grid)
  b <- max(grid)
  p <- length(grid)

  # Step 1: Construct tensor-product basis for grid x grid
  grid_expand <- expand.grid(s = grid, t = grid)
  grid_expand$cov <- as.vector(cov)

  # Spline basis for original grid
  B <- splines2::naturalSpline(grid, derivs = 0, intercept = TRUE, knots = grid[c(-1, -p)], Boundary.knots = c(a, b))  # p x p

  # Use Kronecker product to construct Phi (vectorized outer products)
  Phi <- kronecker(B, B)  # (p^2 x p^2) matrix

  # Step 2: Solve for spline coefficients
  if (any(is.na(cov))) {
    valid_idx <- which(!is.na(grid_expand$cov))
    coef_spline <- qr.solve(Phi[valid_idx, ], grid_expand$cov[valid_idx])
  } else {
    coef_spline <- solve(Phi, grid_expand$cov)
  }

  # Step 3: Interpolation on fine grid
  grid_fine <- seq(a, b, length.out = n_intrpl)
  B_fine <- splines2::naturalSpline(grid_fine, derivs = 0, intercept = TRUE, knots = grid[c(-1, -p)], Boundary.knots = c(a, b))  # (n_fine x p)

  # Fine grid interpolation using Kronecker
  Phi_fine <- kronecker(B_fine, B_fine)  # ((n_intrpl^2) x p^2)

  # Step 4: Predict interpolated covariance
  cov_spline <- matrix(Phi_fine %*% coef_spline, nrow = n_intrpl, ncol = n_intrpl)

  return(cov_spline)

  options(warn=1)

  }
