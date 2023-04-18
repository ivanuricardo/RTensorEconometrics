###### Functions
library(plyr)
library(tensor)

#' Unfold an array along a specified dimension
#'
#' This function performs the "unfold" or "flatten" operation on a given array
#' along a specified dimension.
#'
#' @param X Input array
#' @param n Dimension to unfold
#'
#' @return Unfolded array
#'
#' @examples
#' X <- array(1:24, dim = c(2, 3, 4))
#' unfold(X, 2)
#' returns a matrix with dimensions 3 x 8
#'
#' @import plyr
unfold <- function(X, n) {
  
  # Calculate the combined dimension index and number of dimensions
  idx <- prod(dim(X)[-n])  # combined dimension idx
  ndims <- length(dim(X))
  
  # Create a matrix to store the indices of the flattened array
  mat <- matrix(nrow = idx , ncol = ndims-1)
  
  # Compute indices for each dimension in the matrix
  for (i in 1:(ndims-1)) {
    if (i==1) {
      # Indices for first dimension
      a <- dim(X)[-n][1]
      mat[,1] <- rep((1:a), nrow(mat) / a )
    }
    else {
      # Indices for remaining dimensions
      a <- dim(X)[-n][i]
      a <- unlist(lapply((1:a), rep, times = prod(dim(X)[-n][1:(i-1)])))
      mat[,i] <- rep(a, nrow(mat) / length(a))
    }
  }
  
  # Create a matrix to store the unfolded array
  result <- matrix(nrow = dim(X)[n], ncol = idx)
  
  # Fill the result matrix with values from the original array
  for (i in 1:ncol(result)) {
    # Extract indices from the matrix
    indices <- as.list(mat[i,])
    # Extract values from the original array using the specified indices
    result[,i] <- plyr::take(X, (1:ndims)[-n], indices = indices)
  }
  
  # Return the unfolded array
  result
}

#' Invert a tensor
#'
#' Given a tensor A, this function returns the inverse of A, where the inverse is
#' defined as the inverse of the matricization of A along the first half of its
#' dimensions. 
#'
#' @param A a tensor to be inverted
#' @return the inverse of A
#' @examples
#' A <- array(1:8, dim = c(2, 2, 2))
#' tensor_inverse(A) # returns the inverse of A
tensor_inverse <- function(A) {
  dims <- dim(A)  # get the dimensions of A
  num_dims <- length(dims)  # get the number of dimensions in A
  if (num_dims %% 2 == 1) {  # check if num_dims is odd
    stop("Dimensions must be of even length")
  }
  if (num_dims == 2) {  # check if A is a matrix
    return(solve(A))
  }
  n <- num_dims %/% 2  # integer division to get half the length of dims
  nrows <- prod(dims[1:n])  # calculate number of rows in the resulting mat_A
  matricized_A <- matrix(A, nrow = nrows, byrow = FALSE)  # matricize A
  inverted_A <- solve(matricized_A)  # find inverse of matricized_A
  reshaped_A <- array(inverted_A, dim = dims)  # reshape inverted_A into a tensor
  return(reshaped_A)  # return the inverse of A
}


#' Higher Order Ordinary Least Squares (HOOLS) Estimator 
#'
#' @description Calculates the Ordinary Least Squares (OLS) Estimator for a
#' higher order array. The last dimension of `Y` and `X` is assumed to contain
#' the observations by default. IMPORTANT: Assumes that both response and 
#' predictor tensors are of the same dimensions.
#'
#' @param Y A higher order array with dimensions `d1 x d2 x ... x dn`.
#' @param X A higher order array with dimensions `d1 x d2 x ... x dn`.
#' @param obs_dim_Y An integer indicating the dimension containing the
#' observations in `Y` (default is `n`).
#' @param obs_dim_X An integer indicating the dimension containing the
#' observations in `X` (default is `n`).
#'
#' @return Returns the OLS estimator.
#'
#' @examples
#' HOOLS(Y, X)
#' 
#' @importFrom stats solve
#' @importFrom stats %*%
#' @export

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

