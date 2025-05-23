% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data_doc.R
\docType{data}
\name{simulated_nonstagger_example}
\alias{simulated_nonstagger_example}
\title{A simulated example of non-staggered DiD design}
\format{
It contains two data frames. The first data frame 'simulated_nonstagger_data' has 8200 rows and 5 columns:
\describe{
\item{y}{Outcome variable.}
\item{t}{Time index. There are 41 time periods.}
\item{i}{Unit index. There are 200 units.}
}
The second data frame 'simulated_nonstagger_treatment' has 200 rows and 2 columns:
\describe{
\item{i}{Unit index. It also corresponds to the unit index in data frame 'simulated_nonstagger_data'.}
\item{t0}{Treatment Timing. NA indicates control units.}
\item{x1, x2}{Pre-determined covariates.}
}
}
\description{
A simulated example of non-staggered DiD design
}
