##### Functions

#' CP ALS algorithm
#'
#' Computes the CP Decomposition of a three-dimensional array
#' using the Alternating Least Squares (ALS) method. The algorithm
#' initializes factor matrices A, B, and C from the singular values of X
#' and finds the least-squares solution to each of the subproblems in
#' turn until convergence or a maximum number of iterations is reached.
#'
#' @param X A three-dimensional array.
#' @param R The rank of the CP decomposition.
#' @param tol The convergence tolerance of the algorithm (default 1e-5).
#' @param max_iter The maximum number of iterations of the algorithm (default 1000).
#'
#' @return A list of the factor matrices A, B, and C.
#'
#' @import Matrix
#' @import tensorFun
#' @export

library(Matrix)
library(tensorFun)

cp_als <- function(X, R, tol = 1e-10, max_iter = 1000) {
  # Initialize the three dimensions
  i <- dim(X)[1]
  j <- dim(X)[2]
  k <- dim(X)[3]
  # Initialize from singular values 
  init_mat <- function(X, R, idx) {
    svd_X <- svd(tensorFun::unfold(X, idx))
    sing_X <- diag(svd_X$d) %*% t(svd_X$v)
    return(sing_X[,1:R])
  }
  A_init <- init_mat(X, R, 1)
  B_init <- init_mat(X, R, 2)
  C_init <- init_mat(X, R, 3)
  
  # Initial Estimate and convergence plot
  X_hat <- reconstruct_cp(A_init, B_init, C_init,r = R)
  converge_plot <- sqrt(tensor_ip(X-X_hat))
  
  # Auxiliary Functions
  norm_vec <- function(vec) {norm(as.matrix(vec), type = '2')}
  als_solution <- function(mat_1, mat_2, X, idx) {
    V <- ((t(mat_2) %*% mat_2) * (t(mat_1) %*% mat_1))
    als_hat <- tensorFun::unfold(X, idx) %*% KhatriRao(mat_2, mat_1) %*% solve(V)
    lambdas <- apply(als_hat, 2, norm_vec)
    als_norm <- sweep(als_hat, 2, lambdas, FUN = "/")
    return(list(solution = als_norm, lambdas = lambdas))
  }
  
  A_sol <- als_solution(B_init, C_init, X, 1)
  B_sol <- als_solution(A_sol$solution, C_init, X, 2)
  C_sol <- als_solution(A_sol$solution, B_sol$solution, X, 3)
  
  # Define convergence condition
  convergence_condition <- function(c_list) {
    abs(c_list[length(c_list)] - c_list[length(c_list) - 1]) > tol
  }
  
  # Run iterations until convergence or maximum iterations reached
  iter <- 0
  converged <- FALSE
  
  while((iter < max_iter) && (!converged)) {
    iter <- iter+1
    
    A_sol <- als_solution(B_sol$solution, C_sol$solution, X, 1)
    B_sol <- als_solution(A_sol$solution, C_sol$solution, X, 2)
    C_sol <- als_solution(A_sol$solution, B_sol$solution, X, 3)
    
    X_hat <- reconstruct_cp(A_sol$solution,B_sol$solution,C_sol$solution,r=R, lambda = C_sol$lambdas)
    estimate <- sqrt(tensor_ip(X-X_hat))
    converge_plot[iter+1] <- estimate
    
    if(convergence_condition(converge_plot) == TRUE) {
      converged <- TRUE
    }
  }
  
  X_norm <- sqrt(tensor_ip(X))
  res_norm <- sqrt(tensor_ip(X-X_hat))
  norm_perc <- (1 - (res_norm/X_norm))*100
  
  return(list(A = A_sol, B = B_sol, C = C_sol, est = X_hat, lambdas = C_sol$lambdas, 
              norm_percent = norm_perc, converge_series = converge_plot))
}

#' Reconstruct CP tensor from factor matrices
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

#' tensor_ip
#' Compute the inner product of two three-dimensional tensors.
#'
#' @param X A three-dimensional tensor.
#' @param Y An optional three-dimensional tensor to take the inner
#' product with. If not provided, defaults to \code{X}.
#' @return The inner product of \code{X} and \code{Y}.
#' @details The norm of a three-dimensional tensor is defined as the
#' square root of the sum of squares of all elements. This function
#' computes the inner product of two three-dimensional tensors,
#' analogous to the Frobenius norm for matrices.
#' @examples
#' X <- array(1:27, c(3, 3, 3))
#' tensor_ip(X) # returns 3555
tensor_ip <- function(X, Y = X) {
  return(sum(X * Y))
}

#' Tensor-based Higher Order Singular Value Decomposition (HOSVD)
#'
#' This function performs the Tensor-based Higher Order Singular Value Decomposition
#' (HOSVD) on a given tensor \code{X} using a specified number of ranks \code{num_ranks}.
#'
#' @param X A tensor to be decomposed.
#' @param num_ranks A vector specifying the number of ranks for each mode of \code{X}.
#'
#' @return A list containing:
#' \item{G}{A core tensor.}
#' \item{U}{A list of factor matrices for each mode of \code{X}.}
#' \item{X_hat}{The estimated tensor using the HOSVD.}
#' \item{resid}{The Frobenius norm of the difference between \code{X} and \code{X_hat}.}
#'
#' @import tensorFun
#' @importFrom stats svd
#'
#' @examples
#' X <- array(1:24, dim=c(2,3,4))
#' tucker_hosvd(X, num_ranks=c(2,2,2))
#'
#' @seealso \code{\link{tensor}}, \code{\link{tensor_ip}}
#'
#' @export
tucker_hosvd <- function(X, num_ranks) {
  num_modes <- length(dim(X))
  U_list <- vector("list", num_modes)
  
  # Factor Matrices
  for (m in 1:num_modes) {
    temp_mat <- tensorFun::unfold(X, m)
    U_list[[m]] <- svd(temp_mat, nu = num_ranks[m])$u
  }
  
  # Core Tensor
  G <- X
  for (n in 1:num_modes) {
    G <- tensor(G, t(U_list[[n]]), alongA = 1, alongB = 2)
  }
  
  # Estimate
  X_hat <- G
  for (o in 1:num_modes) {
    X_hat <- tensor(X_hat, U_list[[o]], alongA = 1, alongB = 2)
  }
  
  # Compute residual and return results
  fnorm_X <- sqrt(tensor_ip(X - X_hat))
  return(list(G = G, U = U_list, X_hat = X_hat, resid = fnorm_X))
}








