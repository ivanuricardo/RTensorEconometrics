library(tensor)
library(vars)
library(tsDyn)
library(forecast)
var_ols <- function(A){
  dimensions_A <- dim(A)
  original_A <- A[-1,]
  lagged_A <- A[1:(dimensions_A[1]-1),]
  ols_hat <- solve(t(lagged_A)%*%lagged_A) %*% (t(lagged_A)%*%original_A)
  return(t(ols_hat))
}
# Observations are listed last
HOOLS <- function(Y, X, obs_dim_Y = length(dim(Y)), obs_dim_X = length(dim(X))) {
  
  # Compute numerator and denominator tensors
  numerator <- tensor(Y, X, obs_dim_Y, obs_dim_X)   
  denominator <- tensor(X, X, obs_dim_X, obs_dim_X)
  
  # Invert the denominator tensor
  inverted_den <- tensor_inverse(denominator)
  
  # Compute HOOLS estimator
  number_dims <- length(dim(numerator))
  if (number_dims == 2) {
    # Case when numerator is a matrix
    return(numerator %*% inverted_den)
  } else {
    # Case when numerator is a tensor of higher order
    ols_hat <- tensor(numerator, inverted_den, (number_dims/2 + 1):number_dims,
                      1:(number_dims/2))
    return(ols_hat)
  }
}
# observations are listed first
HOOLS2 <- function(Y, X, obs_dim_Y = 1, obs_dim_X = 1) {
  
  # Compute numerator and denominator tensors
  numerator <- tensor(Y, X, obs_dim_Y, obs_dim_X)
  denominator <- tensor(X, X, obs_dim_X, obs_dim_X)
  
  # Invert the denominator tensor
  inverted_den <- tensor_inverse(denominator)
  
  # Compute HOOLS estimator
  number_dims <- length(dim(numerator))
  if (number_dims == 2) {
    # Case when numerator is a matrix
    return(numerator %*% inverted_den)
  } else {
    # Case when numerator is a tensor of higher order
    ols_hat <- tensor(inverted_den, numerator, (number_dims/2 + 1):number_dims,
                      1:(number_dims/2))
    return(ols_hat)
  }
}

vec.parameters <- matrix(c(0.2, 0.1, -0.2, -0.4, 0.1, 0.2, 0.5, -0.1,
                           0.4, 0.2, 0.6, -0.3, -0.1, 0.3, -0.11, 0.2,
                           0.1, 0.2, 0.5, 0.15, -0.5, -0.1, -0.1, 0.2,
                           0.23, -0.13, -0.1, 0.4, -0.3, 0.1, -0.3,
                           -0.12, -0.3, 0.4, -0.4, 0.4), nrow = 6, ncol = 6)
mat.parameters <- array(vec.parameters, dim = c(3,2,3,2))
vec.data <- VAR.sim(vec.parameters,n=100, include = 'none')

mat.data <- array(vec.data, dim = c(100,3,2))

obs_first_x <- mat.data[2:100,,]

obs_first_y <- mat.data[1:99,,]

var_ols(vec.data)

HOOLS2(obs_first_y, obs_first_x)


########## Simulate Matrix valued data


