# Create test tensor of dimension 3 x 4 x 2
test_tensor <- array(1:24, dim = c(3,4,2))

test_that("unfolding_output", {
  # Flatten along each dimension
  flat1 <- unfold(test_tensor, 1)
  flat2 <- unfold(test_tensor, 2)
  flat3 <- unfold(test_tensor, 3)
  
  # Use output from Kolda & Bader as test comparison
  expected1 <- matrix(1:24, nrow = 3)
  expected2 <- matrix(c(1,4,7,10,2,5,8,11,3,6,9,12,13,16,19,22,14,17,20,23,15,18,21,
                        24), nrow = 4)
  expected3 <- matrix(1:24, nrow = 2, byrow = TRUE)
  
  # Compare actual and expected
  expect_equal(flat1, expected1)
  expect_equal(flat2, expected2)
  expect_equal(flat3, expected3)
})

test_that("unfolding_type", {
  expected <- unfold(test_tensor, 1)
  expect_type(expected,"integer")
})
