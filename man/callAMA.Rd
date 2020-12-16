\name{callAMA}
\alias{callAMA}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Anderson-Moore Algorithm}
\description{
	Calls the Anderson-Moore Algorithm for the given inputs.  callSparseAim is a wrapper to the C and FORTRAN functions which carry out the computation.
}
\usage{
callAMA(cofh, neq, leads, lags, qmax)
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
\details{
Note: for details about what the H, B, Q, S, script A, script B etc. matrices are, please see [Anderson
2010]. You will need to separately load in the data file in the form of
the "H" matrix that is, the
matrix [H_-T ... H_0 ... H_theta]. For details on how to obtain this matrix given a model and
data, see the papers referenced and use the genHmat function. This utility may be added to this
R package in the future. Output will be a list of "list" objects i.e. the B, S, etc. matrices will be
returned as one-dimensional list objects, embedded within a list. It is pretty straightforward to get
these into matrix form (use outObject[k][[1]] to get the kth output where outObject is the result of
callSparseAim). This applies to the B (structural coefficients) and Q (asymptotic linear constraints)
matrices. Alternatively, you can also use the genBmat and genQmat functions to get these matrices
directly, as matrix objects. You can also obtain the S (observable structure) matrix as well as the
stochastic transition matrices, labelled script A and script B in [Anderson 2010].	
}
\value{
	The following list object: [cofh, hrows, hcols, neq, leads, lags, nstate, qmax, ptr ,cofb, qmat]
	Where:
	The first six elements (cofh, ..., lags) are the same as the inputs
	nstate = leads + lags + 1
	qmax is same as on input
	ptr is a pointer to the integer corresponding to the error code; 0 denotes successful output. See getReturnCode() to "translate".
	cofb = reduced form coefficients matrix
	qmat = asymptotic linear constraints matrix
}
\references{
	Gary S. Anderson. "A reliable and computationally efficient algorithm for imposing the saddle point property in dynamic models." Journal of Economic Dynamics & Control, 2010. \cr

	Gary Anderson and George Moore. "An Efficient Procedure for Solving Linear Perfect Foresight Models."  Unpublished manuscript, Board of Governors of the Federal Reserve System. 1983.

}
\author{ Gary Anderson and Aneesh Raghunandan
}
\note{
	It is preferable to ensure that the matrix H is as you want it before passing it as an argument to this function (and to genBmat, genQmat, etc.)
	However, if you wish to simply load the matrix H from a text file with neq*(nlag+nlead+1) lines, make sure to store the entries column-wise as this is R's default method of storage.

	If you are only interested in, for example, the B matrix or the S matrix, then it is not necessary to call this function.  Use genBmat or genScof instead--both call this function as well.
	To extract a matrix from the jth element of the list returned by this function, use result[j][[1]] (where result is the name of the list).
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
	\code{\link{genBmat}}, \code{\link{genScof}}, \code{\link{genQmat}}, \code{\link{genHmat}}, \code{\link{getStochTrans}}, \code{\link{getFactorMatrices}}
}


\examples{
	hmat = read.table(system.file("extdata/Lcofbob2.cof",package="AMA"))
	hrows = 4
	hcols = 12
	neq = 4
	leads = 1
	lags = 1
	output <- callAMA(hmat, neq, leads, lags)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ "models" }
\keyword{ "& algebra"}% __ONLY ONE__ keyword per line
	






%% ~Make other sections like Warning with \section{Warning }{....} ~



% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ "models" }
\keyword{ "& algebra"}% __ONLY ONE__ keyword per line
