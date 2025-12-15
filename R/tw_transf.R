#' Univariate two-way transformation
#' @description The \code{tw_transf} function is used to apply two-way transformation on a single variable.
#'
#' @param x a numeric vector the two-way transformation is supposed to be performed on.
#' @param t a vector of integers that indicates the time.
#' @param i a vector that indicates the unit.
#'
#' @return The \code{tw_transf} function returns two-way transformed data which are reordered by unit and time.
#' @export
#'
#' @examples
#' set.seed(100)
#' x <- rnorm(20)
#' t <- rep(1:5, each = 4)
#' i <- rep(letters[1:4], times=5)
#' df_transf <- tw_transf(x,t,i)
tw_transf <- function(x, t, i) {

  options(warn=-1)

  # check conditions
  if (!is.numeric(x)) stop("The vector 'x' should be numeric.")
  if (!is.numeric(t)) stop("The vector 't' should be numeric.")
  if (!all(table(t,i)==1)) stop("The data should be balanced.")

  # create a data frame
  df <- data.frame(x=x, t=t, i=i)

  # change the data frame to wide format
  mat               <- as.matrix(stats::reshape(df, idvar = "i", timevar = "t", direction = "wide"))
  rownames(mat)     <- NULL
  rownames          <- mat[, 1]  # extract idvar and use it as row names of mat
  mat               <- mat[, -1]
  storage.mode(mat) <- "numeric"

  # apply two-way transformation
  row_means  <- base::rowMeans(mat)
  col_means  <- base::colMeans(mat)
  total_mean <- mean(mat)
  mat_transf <- sweep(mat, 1, row_means, "-") + total_mean
  mat_transf <- sweep(mat_transf, 2, col_means, "-")

  # convert it back to long format
  colnames <- gsub("x.", "", colnames(mat))
  transf   <- data.frame(
    x = as.vector(t(mat_transf)),
    t = rep(as.numeric(colnames), length(rownames)),
    i = rep(rownames,each=length(colnames)))

  # reorder the transformed data frame
  transf <- transf[order(transf$i, transf$t), ]

  # return the transformed data frame
  return(transf)

  options(warn=1)
}
