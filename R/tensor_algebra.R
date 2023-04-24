ttm<-function(tnsr,mat,m=NULL){
	stopifnot(is.matrix(mat))
	if(is.null(m)) stop("m must be specified")
	mat_dims <- dim(mat)
	modes_in <- tnsr@modes
	stopifnot(modes_in[m]==mat_dims[2])
	modes_out <- modes_in
	modes_out[m] <- mat_dims[1]
	tnsr_m <- rs_unfold(tnsr,m=m)@data
	retarr_m <- mat%*%tnsr_m
	rs_fold(retarr_m,m=m,modes=modes_out)
}

#'Tensor Times List
#'
#'Contracted (m-Mode) product between a Tensor of arbitrary number of modes and a list of matrices. The result is folded back into Tensor.
#'@name ttl
#'@rdname ttl
#'@aliases ttl
#'@details Performs \code{ttm} repeated for a single Tensor and a list of matrices on multiple modes. For instance, suppose we want to do multiply a Tensor object \code{tnsr} with three matrices \code{mat1}, \code{mat2}, \code{mat3} on modes 1, 2, and 3. We could do \code{ttm(ttm(ttm(tnsr,mat1,1),mat2,2),3)}, or we could do \code{ttl(tnsr,list(mat1,mat2,mat3),c(1,2,3))}. The order of the matrices in the list should obviously match the order of the modes. This is a common operation for various Tensor decompositions such as CP and Tucker. For the math on the m-Mode Product, see Kolda and Bader (2009).
#'@export
#'@param tnsr Tensor object with K modes
#'@param list_mat a list of matrices
#'@param ms a vector of modes to contract on (order should match the order of \code{list_mat})
#'@return Tensor object with K modes
#'@seealso  \code{\link{ttm}}
#'@note The returned Tensor does not drop any modes equal to 1.
#'@references T. Kolda, B. Bader, "Tensor decomposition and applications". SIAM Applied Mathematics and Applications 2009.
#'@examples
#'tnsr <- new("Tensor",3L,c(3L,4L,5L),data=runif(60))
#'lizt <- list('mat1' = matrix(runif(30),ncol=3),
#' 'mat2' = matrix(runif(40),ncol=4),
#' 'mat3' = matrix(runif(50),ncol=5))
#'ttl(tnsr,lizt,ms=c(1,2,3))
ttl<-function(tnsr,list_mat,ms=NULL){
	if(is.null(ms)||!is.vector(ms)) stop ("m modes must be specified as a vector")
	if(length(ms)!=length(list_mat)) stop("m modes length does not match list_mat length")
	num_mats <- length(list_mat)
	if(length(unique(ms))!=num_mats) warning("consider pre-multiplying matrices for the same m for speed")
	mat_nrows <- vector("list", num_mats)
	mat_ncols <- vector("list", num_mats)
	for(i in 1:num_mats){
	mat <- list_mat[[i]]
	m <- ms[i]
	mat_dims <- dim(mat)
	modes_in <- tnsr@modes
	stopifnot(modes_in[m]==mat_dims[2])
	modes_out <- modes_in
	modes_out[m] <- mat_dims[1]
	tnsr_m <- rs_unfold(tnsr,m=m)@data
	retarr_m <- mat%*%tnsr_m
	tnsr <- rs_fold(retarr_m,m=m,modes=modes_out)
	}
	tnsr
}