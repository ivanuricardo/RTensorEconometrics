% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tensor_regressions.R
\name{MAR}
\alias{MAR}
\title{Matrix AutoRegressive Model (MAR)}
\usage{
MAR(Y, p)
}
\arguments{
\item{Y}{A tensor of response variables.}

\item{p}{Number of lags to take}
}
\value{
A tensor for the MAR estimator
}
\description{
The MAR estimator is a generalization of the ordinary least squares (OLS) estimator to matrix-valued time series.
}
\details{
This function computes the MAR estimator for a matrix-valued regression model. The MAR estimator is defined
as the solution to the equation X = A_1 X_1 + A_2 X_2 + ... + A_p X_r, where X is a tensor of predictor variables,
X_1, X_2, ..., X_r are subtensors of X obtained by fixing the values of certain indices, and A_1, A_2, ..., A_p
are tensors of MAR estimator.
}
\examples{
# Generate random tensor data
library(TensorEconometrics)
Y <- rand_tensor(c(100, 3,4))
X <- rand_tensor(c(100, 3,2))

# Compute MAR estimator with observations along the first dimension
MAR(Y, X, p)

}
