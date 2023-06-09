set.seed(20230609)
X <- rnorm_tnsr(c(6,5,100))
Y <- rnorm_tnsr(c(6,5,100))
R <- c(5,4,3,2)
convThresh <- 1e-04
max_iter <- 400
idx <- 4

xt_regression <- function (X, Y, init_list, idx) {
  omitted_tensor <- ttm(ttm(ttm(init_list[[1]], init_list[[(4-idx)]], (3-idx)), 
                            init_list[[4]], 3), init_list[[5]], 4)
  H <- ttt(X, omitted_tensor, alongA = (4-idx), alongB = (3-idx))
  unfolded_H <- unfold(H, row_idx = c(2,3), col_idx = c(1,4,5))@data
  vec_U <- solve(unfolded_H %*% t(unfolded_H)) %*% unfolded_H %*% vec(Y)
  return(matrix(vec_U, ncol = omitted_tensor@modes[idx]))
}

yt_regression <- function (X, Y, init_list, idx) {
  # This alternates between solving for U3 then solving for U4
  # In order for this to work, we then first omit U3 (include U4) then omit
  # U4 (include U3). Additionally, first multiply over 4 then 3
  # NOTE: init_list includes the core tensor, so shift idx by one
  omitted_tensor <- ttm(ttm(ttm(init_list[[1]], init_list[[2]], 1),
                            init_list[[3]], 2), init_list[[(8-idx)]], (7-idx))
  O <- ttt(X, omitted_tensor, alongA = c(2,3), alongB = c(1,2))
  unfold_dim <- idx-1
  unfolded_O <- unfold(O, row_idx = unfold_dim,
                       col_idx = setdiff(1:3, unfold_dim))@data
  V <- solve(unfolded_O %*% t(unfolded_O)) %*% unfolded_O %*%
    t(unfold(Y, row_idx = unfold_dim, col_idx = setdiff(1:3, unfold_dim))@data)
  return(t(V))
}

core_regression <- function(X, Y, init_list) {
  Xstar <- ttm(ttm(X, t(init_list[[2]]), 2), t(init_list[[3]]), 3)
  Ystar <- ttm(ttm(Y, MASS::ginv(init_list[[4]]), 2),
               MASS::ginv(init_list[[5]]), 3)
  unfolded_Xstar <- unfold(Xstar, row_idx = 1, col_idx = c(2,3))@data
  unfolded_Ystar <- unfold(Ystar, row_idx = 1, col_idx = c(2,3))@data
  G <- solve(crossprod(unfolded_Xstar)) %*% 
    crossprod(unfolded_Xstar, unfolded_Ystar)
  return(G)
}

tuck_conv <- function(pre_init_list, init_list, idx, convThresh) {
  fnorm_conv <- norm(pre_init_list[[idx+1]] - init_list[[idx+1]])
  if(fnorm_conv < convThresh) return(TRUE)
  return(FALSE)
}

tucker_regression <- function(Y, X, R, convThresh = 1e-05, 
                              max_iter = 400, seed = 0) {
  if (seed > 0) set.seed(seed)
  
  # Generate initial random tensor and random CP decomposition
  init_B <- rand_tensor(c(X@modes[-1], Y@modes[-1]))
  # This is a list of a core tensor plus factor matrices
  init_list <- list(
      rand_tensor(R),
      matrix(rnorm(X@modes[2] * R[1]), nrow = X@modes[2]),
      matrix(rnorm(X@modes[3] * R[2]), nrow = X@modes[3]),
      matrix(rnorm(Y@modes[2] * R[3]), nrow = Y@modes[2]),
      matrix(rnorm(Y@modes[3] * R[4]), nrow = Y@modes[3])
  )
  
  converged <- FALSE
  num_iter <- 0
  while (num_iter < max_iter) {
    num_iter <- num_iter + 1
    
    for (idx in 1:init_B@num_modes) {
      if (idx < (init_B@num_modes/2 + 1)) {
        x_results <- xt_regression(init_list = init_list, Y = Y, X = X,
                                   idx = idx)
        pre_init_list <- init_list
        init_list[[idx+1]] <- x_results
        if (tuck_conv(pre_init_list, init_list, idx, convThresh)) {
          converged <- TRUE
          break
        }
      } else {
        y_results <- yt_regression(init_list = init_list, Y = Y, X = X,
                                         idx = idx)
        pre_init_list <- init_list
        init_list[[idx+1]] <- y_results
        if (tuck_conv(pre_init_list, init_list, idx, convThresh)) {
          converged <- TRUE
          break
        }
      }
    }
    flattened_core_update <- core_regression(X = X, Y = Y, init_list = init_list)
    pre_init_list <- init_list
    core_update <- as.tensor(array(flattened_core_update, dim = R))
    init_list[[1]] <- core_update
    
    fnorm_core <- fnorm(pre_init_list[[1]] - init_list[[1]])
    if (fnorm_core < convThresh) converged <- TRUE
    
    if (converged) break  # Exit the loop if converged
  }
  
  return(list(components=init_list, rebuilt_tnsr=tucker_rebuild(init_list)))
}
