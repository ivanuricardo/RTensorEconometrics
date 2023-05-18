test_that("spark_output", {
  # Prepare four matrices, each with spark of 1, 2, 3, 4 respectively
  mat1 <- matrix(c(0,0,0,0,4,1,2,9,-2,-4,1,3), nrow = 4)
  mat2 <- matrix(c(1,0,0,1,-1,0,0,-1))
  mat3 <- matrix(c(1,2,3,2,3,4,3,4,5), nrow = 3)
  mat4 <- diag(3)
  
  # Bonus for completeness
  mat5 <- matrix(c(1,2,3,2,3,4,3,5,7), nrow = 3)
  mat6 <- matrix(c(1,1,0,1,1,0,0,0,1,0,1,1), nrow = 6)
  
  expected1 <- 1
  expected2 <- 2
  expected3 <- 3
  expected4 <- 4
  expected5 <- 3
  expected6 <- 3
    
  result1 <- spark(mat1)
  result2 <- spark(mat2)
  result3 <- spark(mat3)
  result4 <- spark(mat4)
  result5 <- spark(mat5)
  result6 <- spark(mat6)
  
  expect_equal(result1$K, expected1)
  expect_equal(result2$K, expected2)
  expect_equal(result3$K, expected3)
  expect_equal(result4$K, expected4)
  expect_equal(result5$K, expected5)
  expect_equal(result6$K, expected6)
})
