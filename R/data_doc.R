#' A simulated example of staggered DiD design
#' @docType data
#' @name simulated_stagger_example
#'
#' @format It contains two data frames. The first data frame 'simulated_stagger_data' has 4200 rows and 5 columns:
#' \describe{
#'  \item{y}{Outcome variable.}
#'  \item{x1, x2}{Covariates.}
#'  \item{t}{Time index. There are 21 time periods.}
#'  \item{i}{Unit index. There are 200 units.}
#' }
#' The second data frame 'simulated_stagger_treatment' has 200 rows and 2 columns:
#' \describe{
#'  \item{i}{Unit index. It also corresponds to the unit index in data frame 'simulated_stagger_data'.}
#'  \item{t0}{Treatment Timing. It can be -0.5, 0 and 0.5. NA indicates control units.}
#' }
NULL




#' A simulated example of non-staggered DiD design
#' @docType data
#' @name simulated_nonstagger_example
#'
#' @format It contains two data frames. The first data frame 'simulated_nonstagger_data' has 8200 rows and 5 columns:
#' \describe{
#'  \item{y}{Outcome variable.}
#'  \item{x1, x2}{Covariates.}
#'  \item{t}{Time index. There are 41 time periods.}
#'  \item{i}{Unit index. There are 200 units.}
#' }
#' The second data frame 'simulated_nonstagger_treatment' has 200 rows and 2 columns:
#' \describe{
#'  \item{i}{Unit index. It also corresponds to the unit index in data frame 'simulated_nonstagger_data'.}
#'  \item{t0}{Treatment Timing. NA indicates control units.}
#' }
NULL






#' Event study estimates from Lovenheim and Willen (2019)
#' @docType data
#' @name LWdata
#'
#' @format A list, containing 4 objects:
#' \describe{
#'  \item{beta}{Estimates of event study coefficients.}
#'  \item{cov}{Estimates of covariance matrix.}
#'  \item{t0}{Reference time of DiD design.}
#'  \item{paper}{Paper indexing.}
#' }
#' @references
#' Lovenheim, M.F. and Willen, A. (2019). The Long-Run Effects of Teacher Collective Bargaining. American Economic Journal: Economic Policy 11(3), 292– 324
NULL







#' Event study estimates from Bosch and Campos-Vazquez (2014)
#' @docType data
#' @name BCdata
#'
#' @format A list, containing 4 objects:
#' \describe{
#'  \item{beta}{Estimates of event study coefficients.}
#'  \item{cov}{Estimates of covariance matrix.}
#'  \item{t0}{Reference time of DiD design.}
#'  \item{paper}{Paper indexing.}
#' }
#' @references
#' Bosch, M. and Campos-Vazquez, R.M. (2014). The Trade-Offs of Welfare Policies in Labor Markets with Informal Jobs: The Case of the "Seguro Popular" Program in Mexico. American Economic Journal: Economic Policy 6(4), 71–99
NULL








#' Event study estimates from Lafortune et al. (2017)
#' @docType data
#' @name Ldata
#'
#' @format A list, containing 4 objects:
#' \describe{
#'  \item{beta}{Estimates of event study coefficients.}
#'  \item{cov}{Estimates of covariance matrix.}
#'  \item{t0}{Reference time of DiD design.}
#'  \item{paper}{Paper indexing.}
#' }
#' @references
#' Lafortune, J., Rothstein, J. and Schanzenbach, D.W. (2018). School Finance Reform and the Distribution of Student Achievement. American Economic Journal: Applied Economics 10(2), 1–26.
NULL







#' Event study estimates from Gallagher (2014)
#' @docType data
#' @name Gdata
#'
#' @format A list, containing 4 objects:
#' \describe{
#'  \item{beta}{Estimates of event study coefficients.}
#'  \item{cov}{Estimates of covariance matrix.}
#'  \item{t0}{Reference time of DiD design.}
#'  \item{paper}{Paper indexing.}
#' }
#' @references
#' Gallagher, J. (2014). Learning about an Infrequent Event: Evidence from Flood Insurance Take-Up in the United States. American Economic Journal: Applied Economics 6(3), 206–33.
NULL


