##### Functions

#' Reconstruct CP tensor from factor matrices
#' DEPRECIATED
#'
#' Reconstructs the original CP tensor from its factor matrices, i.e., given the factor matrices
#' A, B, and C (corresponding to the three modes of the tensor) and a vector of CP decomposition
#' coefficients lambda, it computes the tensor T as the sum of the outer products of the columns of
#' the factor matrices scaled by the corresponding coefficient.
#'
#' @param A factor matrix of mode 1
#' @param B factor matrix of mode 2
#' @param C factor matrix of mode 3
#' @param r rank of the CP decomposition
#' @param lambda vector of CP decomposition coefficients, defaults to vector of 1's
#' @return Reconstructed CP tensor
reconstruct_cp <- function(A, B, C, r, lambda = rep(1, r)) {
  dim_tens <- c(dim(A)[1], dim(B)[1], dim(C)[1])
  tens <- array(data = 0, dim = dim_tens)
  for(i in 1:r) {
    outer_prod <- lambda[i] * A[,i] %o% B[,i] %o% C[,i]
    tens <- tens + outer_prod
  }
  return(tens)
}

#' Select the rank of a tensor using CP decomposition
#'
#' This function selects the rank of a tensor based on the ratio of the
#' approximation error between different rank values obtained by CP decomposition.
#' 
#' @param tnsr The input tensor.
#' @param max_rank The maximum rank to consider.
#' @return The selected rank.
#' @export
cp_rank_selection <- function(tnsr, max_rank) {
  # Check input types
  if (!inherits(tnsr, "Tensor")) {
    stop("tnsr must be a Tensor object.")
  }
  
  fnorm_fit <- NULL
  
  for (i in 1:max_rank) {
    # Compute CP decomposition with i components
    sim_cp <- cp(tnsr, num_components = i)
    fnorm_fit[i] <- sim_cp$fnorm_resid
  }
  
  # Compute the ratio of approximation errors
  vec_rank_ratio <- fnorm_fit[-max_rank] / fnorm_fit[-1]
  
  # Return the index of the minimum ratio
  return(list(selected_rank = which.min(vec_rank_ratio), rank_vec = vec_rank_ratio))
}

