\name{getReturnCode}
\alias{getReturnCode}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Return Code}
\description{
	Computes the asymptotic linear constraints matrix ("Q") used by the Anderson-Moore Algorithm.  Arguments are identical to those passed to callSparseAim. 
}
\usage{
	getReturnCode(sparseAimObject)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{sparseAimObject}{
	Output object from callSparseAim.
}
}
\details{
	[to be filled in]
}
\value{
	Return code (corresponds to 0 if call to SparseAIM is successful; prints error message otherwise)
}
\references{
	Gary S. Anderson. "A reliable and computationally efficient algorithm for imposing the saddle point property in dynamic models." Journal of Economic Dynamics & Control, 2010.

	Gary Anderson and George Moore. "An Efficient Procedure for Solving Linear Perfect Foresight Models."  Unpublished manuscript, Board of Governors of the Federal Reserve System. 1983.

}

\note{
	Integer corresponding to error code can be extracted directly from the object that callSparseAim returns; i.e. if res is the output of callSparseAim, then use "returnCode = res[9][[1]]".
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
	\code{\link{callAMA}}
}

\examples{
	hmat = read.table(system.file("extdata/Lcofbob2.cof",package="AMA"))
	hrows = 4
	hcols = 12
	neq = 4
	leads = 1
	lags = 1
	output <- callAMA(hmat, neq, leads, lags)
	errorMsg = getReturnCode(output)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ "models" }
\keyword{ "& algebra"}% __ONLY ONE__ keyword per line



