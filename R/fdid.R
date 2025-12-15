#' Estimate the Event Study Coefficients and the Covariance Matrix
#' @description The \code{fdid} function is used to estimate the event study coefficients and the covariance matrix in non-staggered or staggered adoption design.
#' In staggered design, our estimation does not suffer from the so-called negative weighting problem, since we use carefully-chosen non-negative weights to consolidate
#' estimates from each subgroup. See \href{https://arxiv.org/abs/2512.06804}{Fang and Liebl (2025)} for detailed estimation method.
#'
#' @param data a data frame in which the first variable should be the outcome variable, and the latter two variables are time and unit indices. The outcome variable should be numeric.
#' @param treatment a data frame in which the first variable is the unit index; the second variable indicates the reference time, after which the treatment is given, for each unit.
#' NA implies that the unit is never treated. All pre-determined covariates are placed afterwards. Covariates should be numeric.
#'
#' @return The \code{fdid} function returns a list which includes the estimates of event study coefficients and their covariance. In the output, the event time 0 is considered as the
#' reference period. The output is an object of S3 class \code{"fdid"}.
#' @import dplyr
#' @export
#' @references Fang, C. and Liebl, D. (2025). Making Event Study Plots Honest: A Functional Data Approach to Causal Inference. \href{https://arxiv.org/abs/2512.06804}{arXiv:2512.06804}.
#'
#' @seealso \link{tw_transf}, \link{fdid_scb}
#'
#' @examples
#' data(simulated_stagger_example)
#' fdid_est <- fdid(data=simulated_stagger_data, treatment=simulated_stagger_treatment)
#' plot(fdid_est$beta$coef[,2], fdid_est$beta$coef[,1], type="l",
#'      xlab="Event Time", ylab=expression(hat(beta)), family="Times")
fdid <- function(data,
                 treatment)
  {

  options(warn=-1)

  # rename columns of data frames
  colnames(data) <- c("y", "t", "i")
  colnames(treatment)[1:2] <- c("i","t0")

  #check conditions
  if(!is.numeric(data$t)) stop("The time index should be numeric.")
  if(!is.numeric(data$y)) stop("The outcome variable should be numeric.")
  if(!all(na.omit(unique(treatment$t0)) %in% unique(data$t))) stop("The reference time in the data frame 'treatment' should be among the time index in the data frame 'data'.")
  if(!all(na.omit(unique(treatment$t0)) > min(unique(data$t)) & na.omit(unique(treatment$t0)) < max(unique(data$t)))) stop("The reference time in the data frame 'treatment' should be greater/less than the min/max of time index in the data frame 'data'.")
  if(!(NA %in% unique(treatment$t0))) stop("There must be control units with t0=NA as defined in the data frame 'treatment'.")

  # join the two data frames
  data                <- data[order(data$i, data$t), ]
  treatment           <- treatment[order(treatment$i), ]
  treatment[,-(1:2)]  <- apply( treatment[,-(1:2)], 2, function(x) x-mean(x))

  data                <- dplyr::left_join(data, treatment, by="i") # we need to left join them so that we can have a data list with the name of t0
  data_list           <- split(data, ifelse(is.na(data$t0), "NA", data$t0))

  # define a function for estimating beta and cov by each t0
  fdid_nonstagger <- function(t0) {

    ## apply two way transformation on the data frame 'data'
    data_nonstagger        <- rbind(data_list[[as.character(t0)]], data_list[["NA"]])
    data_nonstagger        <- subset(data_nonstagger, select = -c(t0))
    var_names              <- setdiff(colnames(data_nonstagger), c("t","i"))
    t_i_sorted             <- cbind(data_nonstagger$t, data_nonstagger$i)[order(data_nonstagger$i, data_nonstagger$t),]
    t_sorted               <- t_i_sorted[,1]
    i_sorted               <- t_i_sorted[,2]

    data_nonstagger_transf <- data_nonstagger %>% dplyr::mutate(across(y, ~tw_transf(.x,t,i)$x)) %>% dplyr::select(y) %>% dplyr::mutate(t=t_sorted, i=i_sorted)

    ## pre-process the data frame 'treatment'
    treatment_nonstagger          <- treatment[treatment$t0==t0 | is.na(treatment$t0),]
    treatment_nonstagger$d        <- ifelse(is.na(treatment_nonstagger$t0), 0, 1)
    treatment_nonstagger$d_transf <- treatment_nonstagger$d-mean(treatment_nonstagger$d)
    treatment_nonstagger_transf   <- subset(treatment_nonstagger, select=-c(t0,d))
    treatment_nonstagger_transf   <- treatment_nonstagger_transf[order(match(treatment_nonstagger_transf[,"i"], unique(data_nonstagger_transf$i))),]

    ## join the transformed 'data' and pre-processed 'treatment'
    data_nonstagger_transf <- dplyr::left_join(data_nonstagger_transf, treatment_nonstagger_transf, by="i")
    t_unique               <- unique(data_nonstagger_transf$t)

    ## estimate beta
    est_list        <- lapply(t_unique, function(x) solve(crossprod(as.matrix(subset(data_nonstagger_transf[data_nonstagger_transf$t==x,], select=-c(y,t,i)))), crossprod(as.matrix(subset(data_nonstagger_transf[data_nonstagger_transf$t==x,], select=-c(y,t,i))), as.vector(subset(data_nonstagger_transf[data_nonstagger_transf$t==x,], select=y))$y)))
    names(est_list) <- t_unique
    K               <- ncol(data_nonstagger_transf)-3
    gamma           <- cbind(gamma = sapply(est_list, function(x) x[K]), t = t_unique)
    rownames(gamma) <- NULL
    gamma0          <- gamma[gamma[,"t"]==t0, -2]
    beta            <- cbind(beta = gamma[,"gamma"]-gamma0, t = t_unique)

    ## estimate xi_tilde
    if (K != 1){
      xi_tilde           <- cbind(t(sapply(est_list, function(x) x[-K])), t_unique)
      rownames(xi_tilde) <- NULL
      colnames(xi_tilde) <- c(var_names[-1], "t")
    }

    ## estimate covariance of beta
    y_transf_fitted        <- lapply(t_unique, function(x) as.matrix(subset(data_nonstagger_transf[data_nonstagger_transf$t==x,], select=-c(y,t,i))) %*% est_list[[as.character(x)]])
    names(y_transf_fitted) <- t_unique
    resid_transf           <- lapply(t_unique, function(x) y_transf_fitted[[as.character(x)]]-as.vector(subset(data_nonstagger_transf[data_nonstagger_transf$t==x,], select=y))$y)
    names(resid_transf)    <- t_unique

    N_nonstagger     <- length(unique(data_nonstagger_transf$i))
    cov_gamma        <- matrix(unlist(lapply(t_unique, function(x) lapply(t_unique, function(y) solve(1/N_nonstagger*crossprod(treatment_nonstagger_transf$d_transf)) %*% (1/(N_nonstagger-K)* t(treatment_nonstagger_transf$d_transf) %*% diag(x = as.vector(resid_transf[[as.character(x)]])*as.vector(resid_transf[[as.character(y)]]), nrow = N_nonstagger) %*% treatment_nonstagger_transf$d_transf) %*% solve(1/N_nonstagger*crossprod(treatment_nonstagger_transf$d_transf)))))/N_nonstagger, ncol=length(t_unique))
    cov_gamma_gamma0 <- cov_gamma[,which(t_unique==t0)]
    se_gamma0        <- sqrt(cov_gamma[which(t_unique==t0), which(t_unique==t0)])

    cov_beta                       <- (sweep(cov_gamma*N_nonstagger-cov_gamma_gamma0*N_nonstagger, 2, cov_gamma_gamma0*N_nonstagger, FUN="-") + se_gamma0^2*N_nonstagger)/N_nonstagger
    cov_beta[which(t_unique==t0),] <- 0; cov_beta[,which(t_unique==t0)] <- 0
    rownames(cov_beta)             <- colnames(cov_beta) <- t_unique

    ## standard error of xi_tilde
    if (K != 1) {
      est_var_list          <- lapply(t_unique, function(x) (solve(1/N_nonstagger*crossprod(as.matrix(subset(data_nonstagger_transf[data_nonstagger_transf$t==x,], select=-c(y,t,i))))) %*% (1/(N_nonstagger-K)*t(as.matrix(subset(data_nonstagger_transf[data_nonstagger_transf$t==x,], select=-c(y,t,i)))) %*% diag(as.vector(resid_transf[[as.character(x)]]^2), nrow = N_nonstagger) %*% as.matrix(subset(data_nonstagger_transf[data_nonstagger_transf$t==x,], select=-c(y,t,i)))) %*% solve(1/N_nonstagger*crossprod(as.matrix(subset(data_nonstagger_transf[data_nonstagger_transf$t==x,], select=-c(y,t,i))))))/N_nonstagger)
      se_xi_tilde           <- cbind(t(sapply(est_var_list, function(x) sqrt(diag(x)[-K]))), t_unique)
      colnames(se_xi_tilde) <- c(var_names[-1], "t")
    }

    ## return output
    output <- if (K != 1){
      list(beta = list(coef = beta, cov = cov_beta),
           xi_tilde   = list(coef = xi_tilde,   se  = se_xi_tilde))
    } else {
      list(beta = list(coef = beta, cov = cov_beta))
    }
  }

  # run the function above for each t0
  t0_unique                   <- na.omit(unique(treatment$t0))
  fdid_nonstagger_list        <- lapply(t0_unique, fdid_nonstagger)
  names(fdid_nonstagger_list) <- t0_unique

  if (length(t0_unique)==1) {

    # extract the output directly in non-staggered design with only one t0
    fdid_nonstagger_list[[1]]$beta$coef[,"t"]     <- round(fdid_nonstagger_list[[1]]$beta$coef[,"t"]-t0_unique[1], 6)
    colnames(fdid_nonstagger_list[[1]]$beta$coef) <- c("beta", "event_t")
    final_output                                  <- append(fdid_nonstagger_list[[1]], list(t0=0, df=NULL))

  } else {

    # estimate beta in staggered design with multiple t0
    beta_list      <- lapply(1:length(t0_unique), function(x) cbind(beta=fdid_nonstagger_list[[x]]$beta$coef[,"beta"], event_t= round(fdid_nonstagger_list[[x]]$beta$coef[,"t"]-t0_unique[x], 6)))
    t_unique       <- unique(data_list[["NA"]][["t"]])
    event_t_list   <- sapply(beta_list, function(x) x[,"event_t"])
    event_t_unique <- sort(unique(as.vector(event_t_list)))
    N_t0           <- sapply(t0_unique, function(x) sum(treatment$t0==x, na.rm=TRUE) )

    beta_stagger_list <- lapply(beta_list, function(x) cbind(beta = NA, event_t = event_t_unique))
    beta_stagger_list <- Map(function( beta_stagger_mat,beta_mat) {
      match_idx <- match(beta_mat[, "event_t"], beta_stagger_mat[, "event_t"])  # find matching indices
      valid_idx <- !is.na(match_idx)  # remove NA matches
      beta_stagger_mat[match_idx[valid_idx], "beta"] <- beta_mat[valid_idx, "beta"]
      return(beta_stagger_mat)
    }, beta_stagger_list, beta_list)

    beta_stagger <- sapply(beta_stagger_list, function(x) x[,"beta"])
    weights      <- t(apply(t(t(apply(event_t_list,2, function(x) event_t_unique%in%x ))*N_t0),1, function(y) y/sum(y)))
    beta_stagger <- cbind(beta = rowSums(beta_stagger * weights, na.rm=TRUE), event_t=event_t_unique)

    # estimate covariance of beta in staggered design
    cov_beta_list <- lapply(fdid_nonstagger_list, function(x) x$beta$cov)
    cov_beta_list <- lapply(1:length(cov_beta_list), function(x) {
      rownames(cov_beta_list[[x]]) <- colnames(cov_beta_list[[x]]) <- event_t_list[,x]
      return(cov_beta_list[[x]])
    })

    cov_beta_stagger_list <- lapply(cov_beta_list, function(x) matrix(NA, nrow=length(event_t_unique), ncol=length(event_t_unique), dimnames = list(event_t_unique,event_t_unique)))
    replace_matrix_values <- function(A, B) {
      common_rows <- intersect(rownames(A), rownames(B))  # find common row names
      common_cols <- intersect(colnames(A), colnames(B))  # find common column names
      A[common_rows, common_cols] <- B[common_rows, common_cols]  # replace matching values
      return(A)
    }
    cov_beta_stagger_list <- Map(replace_matrix_values, cov_beta_stagger_list, cov_beta_list)

    weights_mat_list <- lapply(1:ncol(weights), function(x) weights[,x] %o% weights[,x])
    cov_beta_stagger <- Reduce(`+`, Map(function(mat1, mat2) mat1 * mat2,
                                        lapply(cov_beta_stagger_list, function(x) replace(x, is.na(x), 0)),
                                        lapply(weights_mat_list, function(x) replace(x, is.na(x), 0))))
    na_indicator     <- Reduce(`&`, lapply(cov_beta_stagger_list, is.na))
    na_indicator     <- ifelse(na_indicator, NA, 1)
    cov_beta_stagger <- cov_beta_stagger*na_indicator

    # collect all estimates of xi_tilde and its standard errors
    if(length(fdid_nonstagger_list[[1]])!=1) {xi_tilde_list <- lapply(fdid_nonstagger_list, function(x) x$xi_tilde)}

    # final output
    final_output <- if(length(fdid_nonstagger_list[[1]])!=1) {
      list(beta=list(coef= beta_stagger, cov=cov_beta_stagger),
           xi_tilde=xi_tilde_list,
           t0=0,
           df=NULL)
      } else {
      list(beta=list(coef= beta_stagger, cov=cov_beta_stagger),
           t0=0,
           df=NULL)
      }
    }

  class(final_output) <- "fdid"
  return(final_output)

  options(warn=1)
  }

