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

- **ttt**: An extension to the ''ttt'' function found in the MATLAB tensor toolbox (Kolda and Bader 2006). Replaces the ''ttl'' function given in rTensor