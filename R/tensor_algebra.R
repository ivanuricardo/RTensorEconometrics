# Tensor Algebra

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
  first_dims <- A@modes[-alongA]
  last_dims <- B@modes[-alongB]
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

#' Calculate the Spark of a Given Matrix
#'
#' The Spark of a matrix is the smallest number of linearly dependent columns. 
#' This function takes a matrix \code{Phi} as input and returns the Spark of 
#' the matrix, \code{K}, and the indices of the columns that make up the 
#' linearly dependent set, \code{columns}.
#'
#' @param Phi a matrix
#' 
#' @return A list with two elements: \code{K}, an integer specifying the 
#'   Spark of the matrix, and \code{columns}, a numeric vector of indices 
#'   specifying the columns that make up the linearly dependent set.
#' 
#' @details This function iteratively chooses subsets of columns of \code{Phi} 
#' and checks if they are linearly dependent. It stops when it finds the 
#' smallest set of linearly dependent columns. If all subsets of columns are 
#' linearly independent, it returns the rank of the matrix plus one and all 
#' columns of the matrix.
#' 
#' @examples
#' Phi <- matrix(c(1, 2, 3, 2, 4, 6, 3, 6, 9), nrow = 3, ncol = 3)
#' spark(Phi)
#' 
#' @export
spark <- function(Phi) {
  # Get the rank of the matrix
  R <- qr(Phi)$rank
  
  # Number of rows and columns of Phi
  N <- ncol(Phi)
  
  # Iterate over the number of columns to select
  for (K in 1:(R + 1)) {
    
    # Check if we have already selected all columns
    if (K > N) {
      break
    }
    
    # Compute the number of possible choices of K columns from N
    numChoices <- choose(N, K)
    
    # Check if the number of choices is reasonable
    if (numChoices > 160000) {
      stop("N choose K too large!")
    }
    
    # Generate all possible choices of K columns out of N
    choices <- combn(1:N, K)
    
    # Iterate over the choices
    for (c in 1:numChoices) {
      choice <- choices[, c]
      Phik <- Phi[, choice]
      
      # Check if the columns are linearly dependent
      if (qr(Phik)$rank < K) {
        return(list(K = K, columns = choice))
      }
    }
  }
  
  # All columns up to rank K are linearly independent
  return(list(K = R + 1, columns = 1:N))
}

#' Random tensor with standard deviation
#'
#' This function generates a random tensor with specified modes and an option
#' for setting the standard deviation.
#'
#' @param modes A vector specifying the dimensions of the tensor.
#' @param sd The standard deviation used to generate random values.
#' @param drop Whether or not modes equal to 1 should be dropped.
#'
#' @return A random tensor with the specified dimensions.
#'
#' @export
rnorm_tnsr <- function(modes = c(3,4,5), sd = 1, drop = FALSE) {
  as.tensor(array(rnorm(prod(modes), sd = sd), dim = modes), drop = drop)
}
