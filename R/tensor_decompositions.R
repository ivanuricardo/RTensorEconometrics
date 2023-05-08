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
#' Additionally, it gives a plot for visual inspection of the fit
#'
#' @param tnsr The input tensor.
#' @param max_rank The maximum rank to consider.
#' @return A list with the selected rank and a vector of the approximation error
#' for each rank.
#' @export
#' @examples
#' # Load data
#' data(traditional_data)
#'
#' # Convert data to tensor
#' traditional_tensor <- as.tensor(traditional_data)
#'
#' # Select the rank of the tensor using CP decomposition
#' cp_rank_selection(traditional_tensor, 20)
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
  
  # Calculate the first and second derivatives of fnorm_fit
  first_deriv <- diff(fnorm_fit)
  second_deriv <- diff(first_deriv)
  
  # Find the index of the maximum value of the second derivative
  kink_index <- which.max(second_deriv)
  
  # Return the rank at the kink
  plot(fnorm_fit)
  return(list(fnorm_fit = fnorm_fit, kink_rank = kink_index + 1))
}


#' Select the rank of each mode of a tensor using Tucker approximation
#'
#' Given a tensor, this function computes the optimal rank for each mode using
#' the SVD of the flattened tensor. WARNING: Can be slow for large tensors.
#' Additionally, for determination of c, see Xia, Xu, and Zhu 2015, Wang et. al 2022.
#'
#' @param tnsr A tensor object
#' @param c Regularization parameter to avoid division by zero. Default is given in 
#' Xia, Xu, and Zhu 2015
#' @return A vector with the estimated rank for each mode
#' @export
tucker_rank_selection <- function(tnsr, c = log(tnsr@modes[1])/(10*tnsr@modes[1])) {
  est_ranks <- NULL
  for (mode in 1:tnsr@num_modes) {
    # Unfold tensor along the current mode
    flattened_tnsr <- unfold(tnsr, mode, setdiff((1:tnsr@num_modes), mode))@data
    
    # Compute the singular values of the flattened tensor
    svd_mode <- svd(flattened_tnsr)$d
    
    # Compute the optimal rank for the current mode
    r_mode <- which.max((svd_mode[1:(tnsr@modes[mode] - 1)] + c) /
                          (svd_mode[2:(tnsr@modes[mode])] + c))
    
    # Store the estimated rank for the current mode
    est_ranks[mode] <- r_mode
  }
  
  # Return the estimated ranks for each mode
  return(est_ranks = est_ranks)
}

