\name{getFactorMatrices}
\alias{getFactorMatrices}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Shock Factor Matrices for Inhomogeneous Solutions}
\description{
	Computes the exogenous shock scaling (phi) and transfer (F) matrices used in the inhomogeneous case (i.e.  
}
\usage{
getFactorMatrices(hmat, bmat, neq, leads, lags)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{hmat}{
	Full H matrix, passed to R columnwise.  The H matrix is the matrix [H_{-T} ... H_0 ... H_{theta}].
	Must have neq rows and neq*(leads+lags+1) columns.
}
  \item{bmat}{
	Full B matrix, as generated either by callSparseAim or genBmat.
}
  \item{neq}{
	Number of equations in the model.  Must equal number of rows in H matrix.
}
  \item{leads}{
	Number of leads in the model.  Must be a positive integer.
}
  \item{lags}{
	Number of lags in the model.  Must be a positive integer.
}
}
\details{
	[to be filled in]
}
\value{
	list object containing phi as first element and F as second element. 
}
\references{
	Gary S. Anderson. "A reliable and computationally efficient algorithm for imposing the saddle point property in dynamic models." Journal of Economic Dynamics & Control, 2010.

	Gary Anderson and George Moore. "An Efficient Procedure for Solving Linear Perfect Foresight Models."  Unpublished manuscript, Board of Governors of the Federal Reserve System. 1983.

}




\seealso{
	\code{\link{callAMA}}, \code{\link{genBmat}}
}


\examples{
	hmat = read.table(system.file("extdata/Lcofbob2.cof",package="AMA"))
	hrows = 4
	hcols = 12
	neq = 4
	leads = 1
	lags = 1
	bmat <- genBmat(hmat, neq, leads, lags)
	hmat = matrix(as.matrix(hmat), neq)
	factorMats <- getFactorMatrices(hmat, bmat, neq, leads, lags)
	phi = factorMats[1][[1]]
	F = factorMats[2][[1]]
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ "models" }
\keyword{ "& algebra"}% __ONLY ONE__ keyword per line


