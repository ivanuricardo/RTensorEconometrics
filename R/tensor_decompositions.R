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
reconstruct_cp <- function(A, B, C, D, r, lambda = rep(1, r)) {
  dim_tens <- c(dim(A)[1], dim(B)[1], dim(C)[1], dim(D)[1])
  tens <- array(data = 0, dim = dim_tens)
  for(i in 1:r) {
    outer_prod <- lambda[i] * A[,i] %o% B[,i] %o% C[,i] %o% D[, i]
    tens <- tens + outer_prod
  }
  return(tens)
}

#'Canonical Polyadic Decomposition
#'
#'Modified version of Canonical Polyadic decomposition. Here, you can
#'specify the type of norm to take for the lambdas, with Euclidean
#'norm being the default
#'@export
#'@details Uses the Alternating Least Squares (ALS) estimation procedure. A progress bar is included to help monitor operations on large tensors.
#'@name cp_modified
#'@rdname cp_modified
#'@aliases cp_modified
#'@param tnsr Tensor with K modes
#'@param num_components the number of rank-1 K-Tensors to use in approximation
#'@param max_iter maximum number of iterations if error stays above \code{tol} 
#'@param tol relative Frobenius norm error tolerance
#'@param norm Type of norm to obtain lambdas. "O" for one norm, "I" 
#'for infinity norm, "F" for frobenius norm, "M" for maximum modulus, 
#'"2" for "spectral" or 2-norm
#'@return a list containing the following \describe{
#'\item{\code{lambdas}}{a vector of normalizing constants, one for each component}
#'\item{\code{U}}{a list of matrices - one for each mode - each matrix with \code{num_components} columns}
#'\item{\code{conv}}{whether or not \code{resid} < \code{tol} by the last iteration}
#'\item{\code{norm_percent}}{the percent of Frobenius norm explained by the approximation}
#'\item{\code{est}}{estimate of \code{tnsr} after compression}
#'\item{\code{fnorm_resid}}{the Frobenius norm of the error \code{fnorm(est-tnsr)}}
#'\item{\code{all_resids}}{vector containing the Frobenius norm of error for all the iterations}
#'}
#'@seealso \code{\link{tucker}}
#'@references T. Kolda, B. Bader, "Tensor decomposition and applications". SIAM Applied Mathematics and Applications 2009.
cp_modified <- function(tnsr, num_components=NULL,max_iter=25, tol=1e-5,
                        norm_type = "2"){
  if(is.null(num_components)) stop("num_components must be specified")
  stopifnot(is(tnsr,"Tensor"))
  if (.is_zero_tensor(tnsr)) stop("Zero tensor detected")
  
  #initialization via truncated hosvd
  num_modes <- tnsr@num_modes
  modes <- tnsr@modes
  U_list <- vector("list",num_modes)
  unfolded_mat <- vector("list",num_modes)
  tnsr_norm <- fnorm(tnsr)
  for(m in 1:num_modes){
    unfolded_mat[[m]] <- rs_unfold(tnsr,m=m)@data
    U_list[[m]] <- matrix(rnorm(modes[m]*num_components), nrow=modes[m], ncol=num_components)
  }
  est <- tnsr
  curr_iter <- 1
  converged <- FALSE
  #set up convergence check
  fnorm_resid <- rep(0, max_iter)
  CHECK_CONV <- function(est){
    curr_resid <- fnorm(est - tnsr)
    fnorm_resid[curr_iter] <<- curr_resid
    if (curr_iter==1) return(FALSE)
    if (abs(curr_resid-fnorm_resid[curr_iter-1])/tnsr_norm < tol) return(TRUE)
    else{ return(FALSE)}
  }	
  #progress bar
  pb <- txtProgressBar(min=0,max=max_iter,style=3)
  #main loop (until convergence or max_iter)
  norm_vec <- function(vec){
    norm(as.matrix(vec), type = norm_type)
  }
  while((curr_iter < max_iter) && (!converged)){
    setTxtProgressBar(pb,curr_iter)
    for(m in 1:num_modes){
      V <- hadamard_list(lapply(U_list[-m],function(x) {t(x)%*%x}))
      V_inv <- solve(V)			
      tmp <- unfolded_mat[[m]]%*%khatri_rao_list(U_list[-m],reverse=TRUE)%*%V_inv
      lambdas <- apply(tmp,2,norm_vec)
      U_list[[m]] <- sweep(tmp,2,lambdas,"/")	
      Z <- .superdiagonal_tensor(num_modes=num_modes,len=num_components,elements=lambdas)
      est <- ttl(Z,U_list,ms=1:num_modes)
    }
    #checks convergence
    if(CHECK_CONV(est)){
      converged <- TRUE
      setTxtProgressBar(pb,max_iter)
    }else{
      curr_iter <- curr_iter + 1
    }
  }
  if(!converged){setTxtProgressBar(pb,max_iter)}
  close(pb)
  #end of main loop
  #put together return list, and returns
  fnorm_resid <- fnorm_resid[fnorm_resid!=0]
  norm_percent<- (1-(tail(fnorm_resid,1)/tnsr_norm))*100
  invisible(list(lambdas=lambdas, U=U_list, conv=converged, est=est, norm_percent=norm_percent, fnorm_resid = tail(fnorm_resid,1),all_resids=fnorm_resid))
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

#' Check for Uniqueness of CP Decomposition
#'
#' Determines whether a CP decomposition is unique based on sufficient and
#' necessary conditions.
#'
#' @param cp_object A CP decomposition object from the 'cptensor' package.
#'
#' @return A list with boolean values for the sufficient and necessary conditions
#'
#' @details The function checks for uniqueness of a CP decomposition based on
#' sufficient and necessary conditions. The sufficient condition checks if the
#' sum of the ranks of each factor matrix is greater than or equal to 2 times
#' the rank of the tensor plus the number of modes minus one. The necessary
#' condition checks if the minimum product of ranks of each factor matrix is
#' greater than or equal to the rank of the tensor.
#'
#' @examples
#' library(rTensor)
#' X <- rand_tensor(c(4,5,6))
#' fit <- cp(X, 3)
#' cp_uniqueness(fit)
#'
#' @importFrom stats qr
#' @export
cp_uniqueness <- function(cp_object) {
  # extract chosen rank and number of modes from original tensor
  R <- length(cp_object$lambdas)
  N <- cp_object$est@num_modes
  
  # Find the sum of k ranks of each factor matrix
  sum_ranks <- sum(sapply(cp_object$U, function(matr) spark(matr)$K))
  suff_cond <- sum_ranks >= 2 * R + (N - 1)
  
  # Find the product of ranks of each factor matrix
  rank_vec <- sapply(cp_object$U, function(matr) qr(matr)$rank)
  prod_rank <- prod(rank_vec) / rank_vec
  necessary_cond <- min(prod_rank) >= R
  
  return(list(suff_cond=suff_cond, necessary_cond=necessary_cond))
}
