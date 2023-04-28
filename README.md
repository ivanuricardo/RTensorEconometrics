# TensorEconometrics
Tensor Methods for Econometrics

This repo can be installed using devtools::install_github("https://github.com/ivanuricardo/TensorEconometrics"). 

Note, this package is meant for use in addition to the rTensor package, and thus
admits the S4 "tensor" class. Note that rTensor has their own implementation of 
CP and Tucker decompositions. 

## Additional Tensor Algebra

- [x] Tensor times Tensor
- [x] Tensor Inverse

## Canonical Polyadic/PARAFAC decomposition

- [x] Decomposition algorithm
- [x] Reconstruction algorithm
- [ ] Higher Order OLS integration

## Tucker Decomposition

- [x] Higher Order SVD
- [ ] Higher Order Orthogonal Iteration
- [ ] Reconstruction algorithm

## Regressions

- [x] Higher Order OLS
- [ ] Matrix Autoregression
- [ ] Tensor Autoregression
- [ ] Tensor-on-Tensor Regression
