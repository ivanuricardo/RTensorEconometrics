
# Only works with tensors of 3 dimensions and observation in the first dim
library(TensorEconometrics)
set.seed(20230517)
Y <- rand_tensor(c(100, 4, 3))
X <- rand_tensor(c(100, 4, 3))
obs_dim_Y <- 1
obs_dim_X <- 1
R <- 6

cp_regression <- function(Y, X, R, obs_dim_X, obs_dim_Y, convThresh = 1e-05, 
                          max_iter = 500, seed = 0) {
  if (seed > 0) set.seed(seed)
  
  # Generate initial random tensor and perform CP decomposition
  init_B <- rand_tensor(c(X@modes[-obs_dim_X], Y@modes[-obs_dim_Y]))
  init_CP <- cp(init_B, R)
  
  # Store the initial CP decomposition factors in a list
  init_list <- init_CP$U
  
  converged <- FALSE
  iter <- 0
  while (!converged || num_iter < max_iter) {
  # Initialize matrices C and D
    Ddims <- c(Y@modes[1] * Y@modes[3], Y@modes[1] * Y@modes[2])
    iter <- iter+1
    for (dims in 1:2) {
      D <- matrix(nrow = Ddims[dims], ncol = 0)
      C <- matrix(nrow = prod(Y@modes), ncol = 0)
      
      for (r in 1:R) {
        # Obtain the CP list without the current dimension
        init_CP_list <- init_list[-dims]
        
        # Extract the rth column from each matrix in the CP list
        rcol_CP <- lapply(init_CP_list, function(m) m[,r])
        
        # Compute the outer product of the extracted vectors
        init_outer <- Reduce(function(x,y) x %o% y, rcol_CP)
        
        # Perform tensor-times-tensor multiplication
        Cr <- ttt(X, as.tensor(init_outer), alongA = (4-dims), alongB = 1)
        
        # Unfold the resulting tensor and append it to matrix C
        unfolded_Cr <- unfold(Cr, c(1,3,4), 2)
        C <- cbind(C, unfolded_Cr@data)
        
        # Obtain the DP list without the current and next dimension
        init_DP_list <- init_list[-(dims+2)]
        
        # Extract the rth column from each matrix in the DP list
        rcol_DP <- lapply(init_DP_list, function(n) n[,r])
        
        # Compute the outer product of the extracted vectors
        init_DP_outer <- Reduce(function(x, y) x %o% y, rcol_DP)
        
        # Perform tensor-tensor-transpose multiplication
        Dr <- ttt(X, as.tensor(init_DP_outer), alongA = 2:3, alongB = 1:2)
        
        # Append the vectorized result to matrix D
        D <- cbind(D, vec(Dr))
      }
      
      # Compute the estimate for B^(1) through ordinary least squares (OLS)
      vec_B1 <- solve(t(C) %*% C) %*% (t(C) %*% vec(Y))
      B1 <- matrix(vec_B1, ncol = R)
      
      # Unfold Y along the dimensions not considered
      Y_unfolded <- unfold(Y, row_idx = setdiff(1:3, (dims+1)),
                           col_idx = (dims+1))@data
      
      # Compute the estimate for B^(3) through OLS
      B3 <- t(solve(t(D) %*% D) %*% t(D) %*% Y_unfolded)
      
      # Update the CP decomposition factors in the list
      init_list[[dims]] <- B1
      init_list[[dims+2]] <- B3
    }
  }
}


library(TensorEconometrics)
set.seed(20230517)
Y <- rand_tensor(c(100, 4, 3))
X <- rand_tensor(c(100, 4, 3))
obs_dim_Y <- 1
obs_dim_X <- 1
R <- 6
max_iter <- 500

cp_regression2 <- function(Y, X, R, obs_dim_X, obs_dim_Y, convThresh = 1e-05, 
                          max_iter = 500, seed = 0) {
  if (seed > 0) set.seed(seed)
  
  # Generate initial random tensor and perform CP decomposition
  init_B <- rand_tensor(c(X@modes[-obs_dim_X], Y@modes[-obs_dim_Y]))
  init_CP <- cp(init_B, R)
  
  # Store the initial CP decomposition factors in a list
  init_list <- init_CP$U
  
  converged <- FALSE
  num_iter <- 0
  while (!converged && num_iter < max_iter) {
  # Initialize matrices C and D
    num_iter <- num_iter + 1
    Ddims <- c(Y@modes[1] * Y@modes[3], Y@modes[1] * Y@modes[2])
    D1 <- matrix(nrow = Ddims[1], ncol = 0)
    D2 <- matrix(nrow = Ddims[2], ncol = 0)
    C1 <- matrix(nrow = prod(Y@modes), ncol = 0)
    C2 <- matrix(nrow = prod(Y@modes), ncol = 0)
    
    # First dimension
    for (r in 1:R) {
      init_CP_list <- init_list[-1]
      rcol_CP <- lapply(init_CP_list, function(m) m[,r])
      init_outer <- Reduce(function(x,y) x %o% y, rcol_CP)
      Cr <- ttt(X, as.tensor(init_outer), alongA = 3, alongB = 1)
      unfolded_Cr <- unfold(Cr, c(1,3,4), 2)
      C1 <- cbind(C1, unfolded_Cr@data)
    }
    vec_B1 <- solve(t(C1) %*% C1) %*% (t(C1) %*% vec(Y))
    B1 <- matrix(vec_B1, ncol = R)
    init_list[[1]] <- B1
    
    # Second Dimension
    for (r in 1:R) {
      init_CP_list <- init_list[-2]
      rcol_CP <- lapply(init_CP_list, function(m) m[,r])
      init_outer <- Reduce(function(x,y) x %o% y, rcol_CP)
      Cr <- ttt(X, as.tensor(init_outer), alongA = 2, alongB = 1)
      unfolded_Cr <- unfold(Cr, c(1,3,4), 2)
      C2 <- cbind(C2, unfolded_Cr@data)
    }
    vec_B2 <- solve(t(C2) %*% C2) %*% (t(C2) %*% vec(Y))
    B2 <- matrix(vec_B2, ncol = R)
    init_list[[2]] <- B2
    
    # Third dimension
    for (r in 1:R) {
      init_DP_list <- init_list[-3]
      rcol_DP <- lapply(init_DP_list, function(n) n[,r])
      init_DP_outer <- Reduce(function(x, y) x %o% y, rcol_DP)
      Dr <- ttt(X, as.tensor(init_DP_outer), alongA = 2:3, alongB = 1:2)
      D1 <- cbind(D1, vec(Dr))
    }
    # Unfold Y along the dimensions not considered
    Y_unfolded <- unfold(Y, row_idx = setdiff(1:3, 2),
                         col_idx = 2)@data
    
    # Compute the estimate for B^(3) through OLS
    B3 <- t(solve(t(D1) %*% D1) %*% t(D1) %*% Y_unfolded)
    init_list[[3]] <- B3
    
    # Fourth dimension
    for (r in 1:R) {
      init_DP_list <- init_list[-4]
      rcol_DP <- lapply(init_DP_list, function(n) n[,r])
      init_DP_outer <- Reduce(function(x, y) x %o% y, rcol_DP)
      Dr <- ttt(X, as.tensor(init_DP_outer), alongA = 2:3, alongB = 1:2)
      D2 <- cbind(D2, vec(Dr))
    }
    # Unfold Y along the dimensions not considered
    Y_unfolded <- unfold(Y, row_idx = setdiff(1:3, 3),
                         col_idx = 3)@data
    
    # Compute the estimate for B^(3) through OLS
    B4 <- t(solve(t(D2) %*% D2) %*% t(D2) %*% Y_unfolded)
    init_list[[4]] <- B4
    
    ## Convergence Condition
    
  }
}