# TensorEconometrics
Tensor Methods for Econometrics

The development version of this package can be installed using devtools::install_github("https://github.com/ivanuricardo/TensorEconometrics").
The stable version can be found via the CRAN repository.

TensorEconometrics was originally built to support my Master Thesis in Economics and Financial Research - Spec. Econometrics.
It has now evolved to be an extension of the rTensor package and provides a multitude of tools for the modern researcher who is interested in tensor/array based data analysis. 
Note that this package is still under development and it is likely to change in the coming months.


TensorEconometrics depends on rTensor, so installation of TensorEconometrics will also install rTensor. 
This means the use of a tensor S4 class applies to all functions provided in this package.
Below you will find the structure of package, including extensions to tensor decompositions, verification techniques, and regression methods incorporating these decompositions.

## Tensor Algebra

- **ttt**: An extension to the ''ttt'' function found in the MATLAB tensor toolbox (Kolda and Bader 2006). Replaces the ''ttl'' function given in rTensor.
- **Spark**: Calculates the Spark (or k-rank) of a given matrix. 
- **rnorm_tnsr**: Generates a random normal tensor with set dimensions and standard deviation.

## Tensor Decompositions

- **cp_modified**: A modification to the CP decomposition given in rTensor. This modification finds lambdas via the frobenius norm (as opposed to the one norm as in rTensor) and removes the progress bar.
- **cp_rank_selection**: Plots the frobenius norm of a cp decomposition for a range of ranks.
- **tucker_rank_selection**: Selects the ranks of each mode for a given tensor based on the eigenvalues of the flattened tensor along the lines of Di Wang, Yao Zheng, Heng Lian & Guodong Li (2022), Lam, Clifford, and Qiwei Yao (2012).
- **cp_uniqueness**: checks uniqueness conditions of CP decomposition based on Kolda and Bader (2006)
- **tucker_rebuild**: Rebuilds a Tucker decomposition from a list of a core tensor and factor matrices.

## Tensor Regressions

- **HOOLS**: Finds the OLS solution for a predictor and response tensor
- **cp_regression**: Finds the OLS solution for a predictor and response tensor with a CP decomposition of the parameters. 
- **tucker_regression**: Finds the OLS solution for a predictor and response tensor with a Tucker decomposition of the parameters.

## Data

- **traditional_data**: Flattened GVAR data
- **tensor_data**: GVAR data in array format