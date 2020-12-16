\name{genScof}
\alias{genScof}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Observable Structure Matrix}
\description{
	Computes the observable structure matrix ("S") using outputs from the Anderson-Moore Algorithm.
}
\usage{
genScof(hmat, bmat, neq, leads, lags)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{hmat}{
	Full H matrix, passed to R columnwise.  The H matrix is the matrix [H_{-T} ... H_0 ... H_{theta}].
}
  \item{bmat}{
	The B (reduced-form coefficients) matrix.
}
  \item{neq}{
	Number of equations in the model.  
}
  \item{leads}{
	Number of leads in the model.  Must be a positive integer.
}
  \item{lags}{
	Number of lags in the model.  Must be a positive integer.
}
}
\value{
	cofS (observable structure matrix)
}
\references{
	Gary S. Anderson. "A reliable and computationally efficient algorithm for imposing the saddle point property in dynamic models." Journal of Economic Dynamics & Control, 2010.

	Gary Anderson and George Moore. "An Efficient Procedure for Solving Linear Perfect Foresight Models."  Unpublished manuscript, Board of Governors of the Federal Reserve System. 1983.

}

\note{
		Alternatively, can be called from the object that callAMA returns; i.e. if res is the output of callAMA, then use "bmat = res[10][[1]]".
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

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
	qmax = 150000
	bmat <- genBmat(hmat, neq, leads, lags)
	Scoef <- genScof(hmat, bmat, neq, leads, lags)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ "models" }
\keyword{ "& algebra"}% __ONLY ONE__ keyword per line


