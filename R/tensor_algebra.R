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