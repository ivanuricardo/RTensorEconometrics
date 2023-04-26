
# Test output is identity
test_that("product_is_identity", {
  library(rTensor)
  set.seed(111111)
  # Two tensor with correct dimensions, one without
  tensor1 <- array(rnorm(36), dim = c(3,2,3,2))
  tensor2 <- array(rnorm(576), dim = c(2,3,4,2,3,4))
  tensor3 <- array(rnorm(24), dim = c(3,4,2))

  inverted1 <- tensor_inverse(tensor1)
  inverted2 <- tensor_inverse(tensor2)
  
  # Find dimension of identity
  dims1 <- dim(tensor1)
  dims2 <- dim(tensor2)
  dim_id1 <- prod(dims1[1:(length(dims1)/2)])
  dim_id2 <- prod(dims2[1:(length(dims2)/2)])
  
  # Create Identities
  identity1 <- diag(dim_id1)
  identity2 <- diag(dim_id2)
  
  # Tensor Product between tensor and inverse
  result1 <- ttt(as.tensor(tensor1), inverted1, alongA = 3:4,
                 alongB = 1:2)
  result2 <- ttt(as.tensor(tensor2), inverted2, alongA = 4:6,
                 alongB = 1:3)
  
  # Test inverse product
  expect_equal(identity1, unfold(result1, c(1,2), c(3,4))@data)
  expect_equal(identity2, unfold(result2, c(1,2,3), c(4,5,6))@data)
  expect_error(tensor_inverse(tensor3))
})
