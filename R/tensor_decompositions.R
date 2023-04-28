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
