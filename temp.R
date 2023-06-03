# Only works with tensors of 3 dimensions and observation in the first dim
devtools::install_github("https://github.com/ivanuricardo/TensorEconometrics")
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
Y <- rand_tensor(c(1000, 4, 3))
X <- rand_tensor(c(1000, 4, 3))
obs_dim_Y <- 1
obs_dim_X <- 1
R <- 6
max_iter <- 500
convThresh <- 1e-05



tmp <- cp_regression2(Y, X, R, obs_dim_X = obs_dim_X, obs_dim_Y = obs_dim_Y)


#### TESTS ####

# In theory, rrr should give the same result as HOOLS if R is large
library(TensorEconometrics)
set.seed(20230518)
X <- rand_tensor(c(1000, 4, 3))
Y <- rand_tensor(c(1000, 4, 3))

obs_dim_X <- 1
obs_dim_Y <- 1
R <- 12

HOOLS_B <- HOOLS(Y, X, 1, 1)
rrr_B <- rrr(X@data, Y@data, R = 12, seed = 20230518)
cp_reg_B <- cp_regression2(Y, X, 12, 1, 1, seed = 20230518)

HOOLS_B@data - cp_reg_B$B
HOOLS_B@data - rrr_B$B

add <- 5