# Tensor Regressions

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
  init_vals <- matrix(stats::rnorm(row_col_prod), nrow = 1)
  
  # Simulate the VAR process
  sim_data <- matrix(0, nrow = row_col_prod, ncol = obs)
  for (i in 2:obs) {
    sim_data[,i] <- parameter_matrix %*% sim_data[,i-1] + stats::rnorm(6)
  }
  sim_data[,1] <- init_vals
  sim_mat <- array(sim_data, dim = c(3,2,obs))
  return(list(matrix_data = sim_mat, vector_data = sim_data))
}

#' Compute the Higher-Order Ordinary Least Squares (HOOLS) estimator
#'
#' The HOOLS estimator is a generalization of the ordinary least squares (OLS) estimator to higher-order tensors.
#' 
#' @param Y A tensor of response variables.
#' @param X A tensor of predictor variables.
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
#' library(TensorEconometrics)
#' Y <- rand_tensor(c(100, 3,4))
#' X <- rand_tensor(c(100, 3,2))
#' 
#' # Compute HOOLS estimator with observations along the first dimension
#' HOOLS(Y, X, obs_dim_Y = 1, obs_dim_X = 1)
#'
#' @importFrom ttt tensor_inverse
#' @export
#' @rdname HOOLS
HOOLS <- function(Y, X) {
  # Check input types
  if (!inherits(Y, "Tensor") || !inherits(X, "Tensor")) {
    stop("Y and X must be tensors")
  }
  full_modes <- c(X@modes[-(X@num_modes)], Y@modes[-(Y@num_modes)])
  flat_Y <- unfold(Y, setdiff(1:Y@num_modes,Y@num_modes), Y@num_modes)@data
  flat_X <- unfold(X, setdiff(1:X@num_modes,X@num_modes), X@num_modes)@data
  
  flat_OLS <- tcrossprod(flat_Y, flat_X) %*% MASS::ginv(tcrossprod(flat_X)) 
  
  return(as.tensor(array(flat_OLS, dim = full_modes)))
}

x1_regression <- function(init_list, X, Y, R) {
  C <- matrix(nrow = prod(Y@modes), ncol = 0)
  omitted_lst <- init_list[-1]
  
  for (r in 1:R) {
    Rvecs <- lapply(omitted_lst, function(m) m[, r])
    outer_vecs <- Reduce(function(x, y) x %o% y, Rvecs)
    Cr <- ttt(X, as.tensor(outer_vecs), alongA = 3,
              alongB = 1)
    unfoldCr <- unfold(Cr, row_idx = 2, col_idx = c(1,3,4))
    C <- cbind(C, t(unfoldCr@data))
  }
  
  vecB1 <- solve(crossprod(C)) %*% crossprod(C, vec(Y))
  unperm_B1 <- matrix(vecB1, nrow = X@modes[2], ncol = R)
  
  # Lambdas
  lambda1 <- apply(unperm_B1, 2, function(x) sqrt(sum(x^2)))
  B1_norm <- sweep(unperm_B1, 2, lambda1, "/")
  
  # Sort matrix in decreasing order
  sorted_idx <- order(lambda1, decreasing = TRUE)
  B1_permuted <- B1_norm[, sorted_idx]
  
  return(list(B1 = as.matrix(B1_permuted), lambdas = lambda1))
}

x2_regression <- function(init_list, X, Y, R) {
  omitted_lst <- init_list[-2]
  C <- matrix(nrow = prod(Y@modes), ncol = 0)
  
  for (r in 1:R) {
    Rvecs <- lapply(omitted_lst, function(m) m[, r])
    outer_vecs <- Reduce(function(x, y) x %o% y, Rvecs)
    Cr <- ttt(X, as.tensor(outer_vecs), alongA = 2, alongB = 1)
    unfoldCr <- unfold(Cr, row_idx = 2, col_idx = c(1,3,4))
    C <- cbind(C, t(unfoldCr@data))
  }
  
  vecB2 <- solve(crossprod(C)) %*% crossprod(C, vec(Y))
  unperm_B2 <- matrix(vecB2, ncol = R)
  
  # Lambdas
  lambda2 <- apply(unperm_B2, 2, function(x) sqrt(sum(x^2)))
  B2_norm <- sweep(unperm_B2, 2, lambda2, "/")
  
  # Sort Matrices
  sorted_idx <- order(lambda2, decreasing = TRUE)
  B2_permuted <- B2_norm[, sorted_idx]
  
  return(list(B2 = as.matrix(B2_permuted), lambdas = lambda2))
}

y1_regression <- function(init_list, X, Y, R) {
  omitted_lst <- init_list[-3]
  D1 <- matrix(nrow = Y@modes[1]*Y@modes[3], ncol = 0)
  
  for (r in 1:R) {
    Rvecs <- lapply(omitted_lst, function(m) m[, r])
    outer_vecs <- Reduce(function(x, y) x %o% y, Rvecs)
    dr <- vec(ttt(X, as.tensor(outer_vecs), alongA = c(2,3), alongB = c(1,2)))
    D1 <- cbind(D1, dr)
  }
  
  unfoldY <- unfold(Y, row_idx = 2, col_idx = c(1,3))
  unperm_B3 <- t(solve(crossprod(D1)) %*% crossprod(D1, t(unfoldY@data)))
  
  # Lambdas
  lambda3 <- apply(unperm_B3, 2, function(x) sqrt(sum(x^2)))
  B3_norm <- sweep(unperm_B3, 2, lambda3, "/")
  
  # Sort Matrices
  sorted_idx <- order(lambda3, decreasing = TRUE)
  B3_permuted <- as.matrix(B3_norm[, sorted_idx])
  
  return(list(B3 = unname(B3_permuted), lambdas = unname(lambda3)))
}

y2_regression <- function(init_list, X, Y, R) {
  omitted_lst <- init_list[-4]
  D2 <- matrix(nrow = Y@modes[1]*Y@modes[2], ncol = 0)
  
  for (r in 1:R) {
    Rvecs <- lapply(omitted_lst, function(m) m[, r])
    outer_vecs <- Reduce(function(x, y) x %o% y, Rvecs)
    dr <- vec(ttt(X, as.tensor(outer_vecs), alongA = c(2,3), alongB = c(1,2)))
    D2 <- cbind(D2, dr)
  }
  
  unfoldY <- unfold(Y, row_idx = 3, col_idx = c(1,2))
  unperm_B4 <- t(solve(crossprod(D2)) %*% crossprod(D2, t(unfoldY@data)))
  
  # Lambdas
  lambda4 <- apply(unperm_B4, 2, function(x) sqrt(sum(x^2)))
  B4_norm <- sweep(unperm_B4, 2, lambda4, "/")
  
  # Sort Matrices
  sorted_idx <- order(lambda4, decreasing = TRUE)
  B4_permuted <- as.matrix(B4_norm[, sorted_idx])
  
  return(list(B4 = unname(B4_permuted), lambdas = unname(lambda4)))
}

init_cp <- function(X, Y, R, obs_dim_X, obs_dim_Y) {
  hools_est <- HOOLS(X=X, Y=Y, obs_dim_Y = 1, obs_dim_X = 1)
  cp_est <- cp_modified(hools_est, num_components = R)
  return(list(cp_est$U[[1]], cp_est$U[[2]], cp_est$U[[3]],
              cp_est$U[[4]]))
}

#' Perform CP regression
#'
#' This function performs CP regression to estimate the tensor B given
#' the input tensor Y and covariate tensor X.
#'
#' @param Y The input tensor.
#' @param X The covariate tensor.
#' @param R The rank of the CP decomposition.
#' @param convThresh The convergence threshold for the CP regression algorithm.
#' @param max_iter The maximum number of iterations.
#' @param seed The random seed for reproducibility.
#'
#' @return A list containing the estimated tensor B, factor matrices, 
#' convergence status, and the number of iterations.
#'
#' @examples
#' Y <- # input tensor
#' X <- # covariate tensor
#' R <- 3
#' obs_dim_X <- # observed dimensions of X
#' obs_dim_Y <- # observed dimensions of Y
#' cp_regression(Y, X, R, obs_dim_X, obs_dim_Y, convThresh = 1e-05,
#'               max_iter = 500, seed = 0)
#'
#' @export
#' @seealso
#' \code{\link{y_regression}}, \code{\link{x_regression}} 
#' \code{\link{conv_cond}}
#' 
#' @export
cp_regression <- function(Y, X, R, convThresh = 1e-04, max_iter = 400,
                          seed = 0) {
  init_list <- init_cp(X=X, Y=Y, R=R, obs_dim_X = obs_dim_X,
                       obs_dim_Y = obs_dim_Y)
  
  converged <- FALSE
  num_iter <- 0
  while (num_iter < max_iter) {
    num_iter <- num_iter + 1
    
    x1_reg <- x1_regression(init_list = init_list, X = X, Y = Y, R = R)
    pre_init_list <- init_list
    init_list[[1]] <- x1_reg$B1
    lambdas <- x1_reg$lambdas
    conv1 <- norm(pre_init_list[[1]] - init_list[[1]], type = "F")
    if (conv1 < convThresh) {
      converged <- TRUE
      break
    }
    
    x2_reg <- x2_regression(init_list = init_list, X = X, Y = Y, R = R)
    pre_init_list <- init_list
    init_list[[2]] <- x2_reg$B2
    lambdas <- x2_reg$lambdas
    conv2 <- norm(pre_init_list[[2]] - init_list[[2]], type = "F")
    if (conv2 < convThresh) {
      converged <- TRUE
      break
    }
    
    y1_reg <- y1_regression(init_list = init_list, X = X, Y = Y, R = R)
    pre_init_list <- init_list
    init_list[[3]] <- y1_reg$B3
    lambdas <- y1_reg$lambdas
    conv3 <- norm(pre_init_list[[3]] - init_list[[3]], type = "F")
    if (conv3 < convThresh) {
      converged <- TRUE
      break
    }
    
    y2_reg <- y2_regression(init_list = init_list, X = X, Y = Y, R = R)
    pre_init_list <- init_list
    init_list[[4]] <- y2_reg$B4
    lambdas <- y2_reg$lambdas
    conv4 <- norm(pre_init_list[[4]] - init_list[[4]], type = "F")
    if (conv4 < convThresh) {
      converged <- TRUE
      break
    }
  }
  
  B <- as.tensor(reconstruct_cp(init_list[[1]], init_list[[2]], init_list[[3]],
                                init_list[[4]], r = R, lambda = lambdas))
  
  return(list(B = B, factor_mat = init_list, converged = converged,
              num_iter = num_iter, lambdas = lambdas))
}

# Tucker Regression

U1_reg <- function (X, Y, init_list) {
  omitted_tensor <- ttm(ttm(ttm(init_list[[1]], init_list[[3]], 2),
                            init_list[[4]], 3), init_list[[5]], 4)
  O <- ttt(omitted_tensor, X, alongA = c(3,4), alongB = c(1,2))
  unfolded_O <- unfold(O, row_idx = 1, col_idx = c(2,3))@data
  V <- unfold(Y, row_idx = 1, col_idx = c(2,3))@data %*% t(unfolded_O) %*%
    MASS::ginv(unfolded_O %*% t(unfolded_O))
  return(V)
}

U2_reg <- function (X, Y, init_list) {
  omitted_tensor <- ttm(ttm(ttm(init_list[[1]], init_list[[2]], 1),
                            init_list[[4]], 3), init_list[[5]], 4)
  O <- ttt(omitted_tensor, X, alongA = c(3,4), alongB = c(1,2))
  unfolded_O <- unfold(O, row_idx = 2, col_idx = c(1,3))@data
  V <- unfold(Y, row_idx = 2, col_idx = c(1,3))@data %*% t(unfolded_O) %*%
    MASS::ginv(unfolded_O %*% t(unfolded_O))
  return(V)
}

U3_reg <- function (X, Y, init_list) {
  # Multiply G x2 U2 x3 U3 x4 U4
  omitted_tensor <- ttm(ttm(ttm(init_list[[1]], init_list[[2]], 1), 
                            init_list[[3]], 2), init_list[[5]], 4)
  H <- ttt(omitted_tensor, X, alongA = 4, alongB = 2)
  unfolded_H <- unfold(H, row_idx = c(3,4), col_idx = c(1,2,5))@data
  vec_U <- vec(Y) %*% t(unfolded_H) %*% MASS::ginv(unfolded_H %*% t(unfolded_H))
  return(matrix(vec_U, nrow = X@modes[1]))
}

U4_reg <- function (X, Y, init_list) {
  # Multiply G x2 U2 x3 U3 x4 U4
  omitted_tensor <- ttm(ttm(ttm(init_list[[1]], init_list[[2]], 1), 
                            init_list[[3]], 2), init_list[[4]], 3)
  H <- ttt(omitted_tensor, X, alongA = 3, alongB = 1)
  unfolded_H <- unfold(H, row_idx = c(3,4), col_idx = c(1,2,5))@data
  vec_U <- vec(Y) %*% t(unfolded_H) %*% MASS::ginv(unfolded_H %*% t(unfolded_H))
  return(matrix(vec_U, nrow = X@modes[2]))
}

core_regression <- function(X, Y, R, init_list) {
  Xstar <- ttm(ttm(X, t(init_list[[2]]), 1), t(init_list[[3]]), 2)
  Ystar <- ttm(ttm(Y, MASS::ginv(init_list[[4]]), 1),
               MASS::ginv(init_list[[5]]), 2)
  G <- HOOLS(Ystar, Xstar)
  return(G)
}

init_est <- function(Y, X, R) {
  hools_est <- HOOLS(Y, X)
  hosvd_est <- hosvd(hools_est, ranks = R)
  return(list(hosvd_est$Z, hosvd_est$U[[1]], hosvd_est$U[[2]], hosvd_est$U[[3]],
              hosvd_est$U[[4]]))
}

#' Perform Tucker regression
#'
#' This function performs Tucker regression on a given tensor, decomposing it
#' into a low-rank tensor and factor matrices.
#'
#' @param Y The target tensor.
#' @param X The predictor tensor.
#' @param R A vector specifying the desired rank of each mode in the Tucker
#'   decomposition.
#' @param convThresh The convergence threshold for the Tucker decomposition.
#' @param max_iter The maximum number of iterations for convergence.
#' @param seed An optional seed for reproducible results.
#'
#' @return A list containing the decomposed components and the rebuilt tensor.
#'
#' @export
tucker_regression <- function(Y, X, R, convThresh = 1e-04, max_iter = 400,
                              seed = 0) {
  if (seed > 0) set.seed(seed)
  
  # Initialize via HOSVD
  init_list <- init_est(X = X, Y = Y, R = R)
  init_B <- tucker_rebuild(init_list)
  
  converged <- FALSE
  num_iter <- 0
  
  while (num_iter < max_iter) {
    num_iter <- num_iter + 1
    
    # Update U1, U2, U3, U4, and G
    U1 <- U1_reg(X, Y, init_list = init_list)
    pre_init_list <- init_list
    init_list[[2]] <- U1
    conv1 <- norm(pre_init_list[[2]] - init_list[[2]], type = "F")
    if (conv1 < convThresh) {
      converged <- TRUE
      break
    }
    U2 <- U2_reg(X, Y, init_list = init_list)
    pre_init_list <- init_list
    init_list[[3]] <- U2
    conv2 <- norm(pre_init_list[[3]] - init_list[[3]], type = "F")
    if (conv1 < convThresh) {
      converged <- TRUE
      break
    }
    
    U3 <- U3_reg(X, Y, init_list = init_list)
    pre_init_list <- init_list
    init_list[[4]] <- U3
    conv3 <- norm(pre_init_list[[4]] - init_list[[4]], type = "F")
    if (conv1 < convThresh) {
      converged <- TRUE
      break
    }
    
    U4 <- U4_reg(X, Y, init_list = init_list)
    pre_init_list <- init_list
    init_list[[5]] <- U4
    conv4 <- norm(pre_init_list[[5]] - init_list[[5]], type = "F")
    if (conv1 < convThresh) {
      converged <- TRUE
      break
    }
    
    G <- core_regression(X, Y, R, init_list = init_list)
    pre_init_list <- init_list
    init_list[[1]] <- as.tensor(G)
    conv5 <- fnorm(pre_init_list[[1]] - init_list[[1]])
    if (conv1 < convThresh) {
      converged <- TRUE
      break
    }
  }
  
  # HOSVD to make factors orthogonal
  rebuild_A <- tucker_rebuild(init_list)
  hosvd_A <- hosvd(rebuild_A, ranks = R)
  mse <- Y - ttt(X, hosvd_A$est, alongA = 2:3, alongB = 1:2)
  
  return(list(G = hosvd_A$Z, U = hosvd_A$U, A = hosvd_A$est, 
              num_iter = num_iter, converged = converged, mse = mse))
}
