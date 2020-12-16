\name{getStochTrans}
\alias{getStochTrans}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Stochastic Transition Matrices}
\description{
	Computes the stochastic transition matrices ("script A" and "script B") used by the Anderson-Moore Algorithm.  Matrices H and S must already be matrix objects of proper dimension.
}
\usage{
getStochTrans(hmat, scof)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{hmat}{
	Full H matrix object.  The H matrix is the matrix [H_{-T} ... H_0 ... H_{theta}]. (i.e. T lags and theta leads)
}
  \item{scof}{
	Matrix S ("observable structure") .
}

}
\details{
	A and B are the solution to [x_{t-T+1} ... x_{t}]' = A[x_{t-T} ...x{t-1}]' + B[e_t + psi*(E[z_{t}|I_{t}] - E[z_{t}|I_{t-1})]}
\value{
	A, B (returned as a "list" object; the first element is A and the second is B.)
}
\references{
	Gary S. Anderson. "A reliable and computationally efficient algorithm for imposing the saddle point property in dynamic models." Journal of Economic Dynamics & Control, 2010. \cr

	Gary Anderson and George Moore. "An Efficient Procedure for Solving Linear Perfect Foresight Models."  Unpublished manuscript, Board of Governors of the Federal Reserve System. 1983.

}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
	\code{\link{callAMA}}, \code{\link{genScof}}, \code{\link{genBmat}}
}


\examples{
	hmat = read.table(system.file("extdata/Lcofbob2.cof",package="AMA"))
	hrows = 4
	hcols = 12
	neq = 4
	hmatrix = matrix((as.matrix(hmat)), neq)
	leads = 1
	lags = 1
	bmat <- genBmat(hmat, neq, leads, lags)
	Scoef <- genScof(hmat, bmat, neq, leads, lags)
	stochMatrices <- getStochTrans(hmatrix, Scoef)
	scriptA = stochMatrices[1][[1]]
	scriptB = stochMatrices[2][[1]]
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ "models" }
\keyword{ "& algebra"}% __ONLY ONE__ keyword per line



