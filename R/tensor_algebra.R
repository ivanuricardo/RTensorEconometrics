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
#' A <- rnorm_tnsr(c(4,5,6))
#' B <- rnorm_tnsr(c(6,8,9))
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

#' Tensor times Vector
#' 
#' Given a tensor and a (column) vector, compute the product. The integer N 
#' specifies the dimension in A along which V is multiplied.
#' 
#' @param A A tensor
#' @param V A column vector
#' @param N Dimension over which to multiply
#' 
#' @return A contracted tensor
#' 
#' @export
ttv <- function(A, V, N) {
  if (is.array(A) || is.vector(A)) A <- as.tensor(A)
  return(ttt(A, as.tensor(V), N, 1))
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
