#' Compute the tensor product of two tensors along specified modes.
#'
#' @param A First input tensor.
#' @param B Second input tensor.
#' @param alongA An integer indicating the mode of A along which to perform
#'   the tensor product. Default is NA.
#' @param alongB An integer indicating the mode of B along which to perform
#'   the tensor product. If NA, set to alongA. Default is NA.
#' 
#' @return A new tensor resulting from the tensor product of A and B along the
#'   specified modes.
#'
#' @examples
#' A <- array(1:24, dim = c(3, 4, 2))
#' B <- array(1:24, dim = c(2, 3, 4))
#' ttt(A, B, alongA = 3, alongB = 1)
#' 
#' @export
ttt <- function(A, B, alongA = NA, alongB = NA) {
  stopifnot(is(A, "Tensor"))
  stopifnot(is(B, "Tensor"))
  
  # Set alongB to alongA if alongB is NA
  alongB <- ifelse(is.na(alongB), alongA, alongB)
  
  # Get dimensions of resulting tensor
  first_dims <- A@modes[-match(A@modes[alongA], A@modes)]
  last_dims <- B@modes[-match(B@modes[alongB], B@modes)]
  full_dims <- c(first_dims, last_dims)
  
  # Compute the tensor product
  amatrix <- unfold(A, setdiff(1:A@num_modes, alongA), alongA)
  bmatrix <- unfold(B, alongB, setdiff(1:B@num_modes, alongB))
  cmatrix <- amatrix@data %*% bmatrix@data
  
  # Return the result as a tensor
  as.tensor(array(cmatrix, dim = full_dims))
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
  if (is.array(A) || is.vector(A)) A <- as.tensor(A)
  if (A@num_modes %% 2 == 1) {
    stop("Dimensions must be of even length")
  }
  if(A@num_modes == 2) {
    return(solve(A@data))
  }
  fin_nrows <- prod(A@modes[1:(length(A@modes)/2)])
  matricized_A <- unfold(A, row_idx = 1:(A@num_modes/2),
                         col_idx = (A@num_modes/2 + 1):A@num_modes)
  inverse_A <- solve(matricized_A@data)
  reshaped_A <- array(inverse_A, dim = A@modes)
  return(as.tensor(reshaped_A))
}

#' Simulate data from a VAR process and estimate the OLS coefficients
#'
#' This function simulates data from a VAR process with coefficients given by a parameter matrix and then
#' estimates the OLS coefficients using matrix algebra.
#' 
#' @param obs An integer specifying the number of observations to simulate.
#' @param parameter_matrix A num_rows * num_cols by num_rows * num_cols parameter matrix for the VAR process.
#' @param num_rows An integer specifying the number of rows in the parameter matrix.
#' @param num_cols An integer specifying the number of columns in the parameter matrix.
#'
#' @return A list containing the simulated data as a 3-mode tensor with dimensions obs x num_rows x num_cols and
#' as a matrix with dimensions obs x (num_rows * num_cols).
#' 
#' @details This function simulates data from a VAR process of the form X_t = A_1 X_{t-1} + ... + A_p X_{t-p} + u_t,
#' where X_t is a num_rows x num_cols matrix of observed variables, A_1, ..., A_p are num_rows x num_cols parameter
#' matrices, p is the order of the VAR process, and u_t is a num_rows x num_cols matrix of errors with independent,
#' identically distributed normal entries. The parameter matrix is given by vec(A_1), vec(A_2), ..., vec(A_p),
#' where vec(A_i) denotes the vectorization of the matrix A_i. The OLS coefficients for the VAR process are then
#' estimated using matrix algebra.
#' 
#' @examples
#' # Simulate data from a VAR process with 2 variables and order 1
#' parameter_matrix <- matrix(c(0.5, 0.3, 0.2, 0.4), nrow = 2)
#' matrix_ols_sim(obs = 100, parameter_matrix = parameter_matrix, num_rows = 2, num_cols = 1)
#'
#' @importFrom stats rnorm
#' @export
#' @rdname matrix_ols_sim
matrix_ols_sim <- function(obs, parameter_matrix, num_rows, num_cols) {
  mat_parameters <- aperm(array(parameter_matrix, dim = c(num_rows,num_cols,num_rows,num_cols)), c(1,3,2,4))
  row_col_prod <- num_rows*num_cols
  
  # Generate random normal variables for the initial values
  init_vals <- matrix(rnorm(row_col_prod), nrow = 1)
  
  # Simulate the VAR process
  sim_data <- matrix(0, nrow = obs, ncol = row_col_prod)
  for (i in 2:obs) {
    sim_data[i,] <- parameter_matrix %*% sim_data[i-1,] + rnorm(6)
  }
  sim_data[1, ] <- init_vals
  sim_mat <- array(sim_data, dim = c(obs,3,2))
  return(list(matrix_data = sim_mat, vector_data = sim_data))
}

#' Compute the Higher-Order Ordinary Least Squares (HOOLS) estimator
#'
#' The HOOLS estimator is a generalization of the ordinary least squares (OLS) estimator to higher-order tensors.
#' 
#' @param Y A tensor of response variables.
#' @param X A tensor of predictor variables.
#' @param obs_dim_Y An integer indicating the dimension of \code{Y} along which observations are stacked.
#' @param obs_dim_X An integer indicating the dimension of \code{X} along which observations are stacked.
#'
#' @return A tensor of HOOLS estimators.
#' 
#' @details This function computes the HOOLS estimator for a tensor regression model. The HOOLS estimator is defined
#' as the solution to the equation X = X_1 B_1 + X_2 B_2 + ... + X_r B_r, where X is a tensor of predictor variables,
#' X_1, X_2, ..., X_r are subtensors of X obtained by fixing the values of certain indices, and B_1, B_2, ..., B_r
#' are tensors of HOOLS estimators.
#' 
#' @examples
#' # Generate random tensor data
#' library(rTensor)
#' Y <- rand_tensor(c(100, 3,4))
#' X <- rand_tensor(c(100, 3,2))
#' 
#' # Compute HOOLS estimator with observations along the first dimension
#' HOOLS(Y, X, obs_dim_Y = 1, obs_dim_X = 1)
#'
#' @importFrom ttt tensor_inverse
#' @export
#' @rdname HOOLS
HOOLS <- function(Y, X, obs_dim_Y = length(dim(Y)), obs_dim_X = length(dim(X))) {
  # Check input types
  if (!inherits(Y, "Tensor") || !inherits(X, "Tensor")) {
    stop("Y and X must be tensors")
  }
  
  # Compute numerator and denominator tensors
  numerator <- ttt(X, Y, obs_dim_Y, obs_dim_X)
  denominator <- ttt(X, X, obs_dim_X, obs_dim_X)
  
  # Invert the denominator tensor
  inverted_den <- tensor_inverse(denominator)
  
  # Compute HOOLS estimator
  number_dims <- length(dim(numerator))
  ols_hat <- ttt(inverted_den, numerator, (number_dims/2 + 1):number_dims, 1:(number_dims/2))
  
  return(ols_hat)
}



