\name{reconstruct_cp}
\alias{reconstruct_cp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Reconstruct CP tensor from factor matrices
}
\description{
#' Reconstructs the original CP tensor from its factor matrices, i.e., given the factor matrices
#' A, B, and C (corresponding to the three modes of the tensor) and a vector of CP decomposition
#' coefficients lambda, it computes the tensor T as the sum of the outer products of the columns of
#' the factor matrices scaled by the corresponding coefficient
}
\usage{
reconstruct_cp(A, B, C, r, lambda = rep(1, r))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{A}{
%%     ~~Describe \code{A} here~~
}
  \item{B}{
%%     ~~Describe \code{B} here~~
}
  \item{C}{
%%     ~~Describe \code{C} here~~
}
  \item{r}{
%%     ~~Describe \code{r} here~~
}
  \item{lambda}{
%%     ~~Describe \code{lambda} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (A, B, C, r, lambda = rep(1, r)) 
{
    dim_tens <- c(dim(A)[1], dim(B)[1], dim(C)[1])
    tens <- array(data = 0, dim = dim_tens)
    for (i in 1:r) {
        outer_prod <- lambda[i] * A[, i] \%o\% B[, i] \%o\% C[, i]
        tens <- tens + outer_prod
    }
    return(tens)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
