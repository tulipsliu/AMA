\name{genBmat}
\alias{genBmat}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Reduced-Form Coefficients Matrix}
\description{
	Computes the reduced-form coefficients matrix ("B") output by the Anderson-Moore Algorithm.  Arguments are identical to those passed to callSparseAim. 
}
\usage{
genBmat(cofh, neq, leads, lags, qmax)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{cofh}{
	Full H matrix, passed to R columnwise.  The H matrix is the matrix [H_{-T} ... H_0 ... H_{theta}].
	Must have neq rows and neq*(leads+lags+1) columns.
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
  \item{qmax}{
	Number of elements to allocate for sparse matrix storage.  Optional input.  Default is neq * leads * neq * (leads + lags).
}
}

\value{
	cofb ("B" matrix)
}
\references{
	Gary S. Anderson. "A reliable and computationally efficient algorithm for imposing the saddle point property in dynamic models." Journal of Economic Dynamics & Control, 2010. \cr

	Gary Anderson and George Moore. "An Efficient Procedure for Solving Linear Perfect Foresight Models."  Unpublished manuscript, Board of Governors of the Federal Reserve System. 1983.

}

\note{
		Alternatively, can be called from the object that callAMA returns; i.e. if res is the output of callAMA, then use "bmat = res[10][[1]]".
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
	\code{\link{callAMA}}
}


\examples{
	hmat <- read.table(system.file("extdata/Lcofbob2.cof",package="AMA"))
	bmat <- genBmat(hmat, 4, 1, 1)


hmat<-genHmat(system.file("extdata/example7.mod",package="AMA"),system.file("extdata/example7params.prm",package="AMA"))
bmat<-genBmat(hmat,4,1,1)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ "models" }
\keyword{ "& algebra"}% __ONLY ONE__ keyword per line

