\name{cp_als}
\alias{cp_als}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
CP ALS algorithm
}
\description{
#' Computes the CP Decomposition of a three-dimensional array
#' using the Alternating Least Squares (ALS) method. The algorithm
#' initializes factor matrices A, B, and C from the singular values of X
#' and finds the least-squares solution to each of the subproblems in
#' turn until convergence or a maximum number of iterations is reached.
}
\usage{
cp_als(X, R, tol = 1e-10, max_iter = 1000)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
%%     ~~Describe \code{X} here~~
}
  \item{R}{
%%     ~~Describe \code{R} here~~
}
  \item{tol}{
%%     ~~Describe \code{tol} here~~
}
  \item{max_iter}{
%%     ~~Describe \code{max_iter} here~~
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
function (X, R, tol = 1e-10, max_iter = 1000) 
{
    i <- dim(X)[1]
    j <- dim(X)[2]
    k <- dim(X)[3]
    init_mat <- function(X, R, idx) {
        svd_X <- svd(tensorFun::unfold(X, idx))
        sing_X <- diag(svd_X$d) \%*\% t(svd_X$v)
        return(sing_X[, 1:R])
    }
    A_init <- init_mat(X, R, 1)
    B_init <- init_mat(X, R, 2)
    C_init <- init_mat(X, R, 3)
    X_hat <- reconstruct_cp(A_init, B_init, C_init, r = R)
    converge_plot <- sqrt(tensor_ip(X - X_hat))
    norm_vec <- function(vec) {
        norm(as.matrix(vec), type = "2")
    }
    als_solution <- function(mat_1, mat_2, X, idx) {
        V <- ((t(mat_2) \%*\% mat_2) * (t(mat_1) \%*\% mat_1))
        als_hat <- tensorFun::unfold(X, idx) \%*\% KhatriRao(mat_2, 
            mat_1) \%*\% solve(V)
        lambdas <- apply(als_hat, 2, norm_vec)
        als_norm <- sweep(als_hat, 2, lambdas, FUN = "/")
        return(list(solution = als_norm, lambdas = lambdas))
    }
    A_sol <- als_solution(B_init, C_init, X, 1)
    B_sol <- als_solution(A_sol$solution, C_init, X, 2)
    C_sol <- als_solution(A_sol$solution, B_sol$solution, X, 
        3)
    convergence_condition <- function(c_list) {
        abs(c_list[length(c_list)] - c_list[length(c_list) - 
            1]) > tol
    }
    iter <- 0
    converged <- FALSE
    while ((iter < max_iter) && (!converged)) {
        iter <- iter + 1
        A_sol <- als_solution(B_sol$solution, C_sol$solution, 
            X, 1)
        B_sol <- als_solution(A_sol$solution, C_sol$solution, 
            X, 2)
        C_sol <- als_solution(A_sol$solution, B_sol$solution, 
            X, 3)
        X_hat <- reconstruct_cp(A_sol$solution, B_sol$solution, 
            C_sol$solution, r = R, lambda = C_sol$lambdas)
        estimate <- sqrt(tensor_ip(X - X_hat))
        converge_plot[iter + 1] <- estimate
        if (convergence_condition(converge_plot) == TRUE) {
            converged <- TRUE
        }
    }
    X_norm <- sqrt(tensor_ip(X))
    res_norm <- sqrt(tensor_ip(X - X_hat))
    norm_perc <- (1 - (res_norm/X_norm)) * 100
    return(list(A = A_sol, B = B_sol, C = C_sol, est = X_hat, 
        lambdas = C_sol$lambdas, norm_percent = norm_perc, converge_series = converge_plot))
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
