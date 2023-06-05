test_that("Tucker Rebuilds", {
  set.seed(20230604)
  
  # Create original tensor
  original_tensor <- rand_tensor(c(6,5,4,3))
  
  # Specify modes
  R <- c(5,4,3,2)
  
  # Tucker decompose
  tucker_tensor <- tucker(original_tensor, R)
  
  tucker_list <- list(tucker_tensor$Z, tucker_tensor$U[[1]],
                      tucker_tensor$U[[2]], tucker_tensor$U[[3]], 
                      tucker_tensor$U[[4]])
  rebuilt_tensor <- tucker_rebuild(tucker_list)
  
  expect_equal(rebuilt_tensor@num_modes, original_tensor@num_modes)
  expect_equal(rebuilt_tensor@modes, original_tensor@modes)
  expect_equal(rebuilt_tensor@data, original_tensor@data, tolerance = 1)
  
  
})