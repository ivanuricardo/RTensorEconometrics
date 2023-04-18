\name{tucker_hosvd}
\alias{tucker_hosvd}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Tucker HOSVD
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
tucker_hosvd(X, num_ranks)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
%%     ~~Describe \code{X} here~~
}
  \item{num_ranks}{
%%     ~~Describe \code{num_ranks} here~~
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
function (X, num_ranks) 
{
    num_modes <- length(dim(X))
    U_list <- vector("list", num_modes)
    for (m in 1:num_modes) {
        temp_mat <- tensorFun::unfold(X, m)
        U_list[[m]] <- svd(temp_mat, nu = num_ranks[m])$u
    }
    G <- X
    for (n in 1:num_modes) {
        G <- tensor(G, t(U_list[[n]]), alongA = 1, alongB = 2)
    }
    X_hat <- G
    for (o in 1:num_modes) {
        X_hat <- tensor(X_hat, U_list[[o]], alongA = 1, alongB = 2)
    }
    fnorm_X <- sqrt(tensor_ip(X - X_hat))
    return(list(G = G, U = U_list, X_hat = X_hat, resid = fnorm_X))
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
