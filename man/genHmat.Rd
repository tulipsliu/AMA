\name{genHmat}
\alias{genHmat}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Structural Coefficients Matrix Scripts}
\description{
	Creates the structural coefficients matrix H, given the model file and parameter values as inputs.
}
\usage{
genHmat(modelFileNameFull, params, wantParamVec = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{modelFileNameFull}{
	String object containing the name of the model file.  Model file must be in MODELEZ syntax.
}
  \item{params}{
	Either:
	-A string object containing the name of the parameters file (One parameter per line, with each line in the format "param = <value>")
	OR:
	-A vector containing the parameters.  Vector must have names associated with each entry (i.e. names(params) should be correct and complete)
  }
  \item{wantParamVec}{
	If TRUE, then the output of the function is a list object containing as its first element the H matrix and as its second element the vector of parameters.  Default is FALSE.
	You should probably set this to TRUE if you're reading in from a parameter file, and leave it as FALSE if you already have a parameter vector.
  }
}
\details{
	Model input file must be in MODELEZ syntax.
}
\value{
	Structural coefficients matrix H.
}
\references{
For the details of MODELEZ syntax, see for example \url{http://www.federalreserve.gov/pubs/oss/oss4/papers/habitMatlab/habitMatlab.html}

For details about the algorithm:

	Gary S. Anderson. "A reliable and computationally efficient algorithm for imposing the saddle point property in dynamic models." Journal of Economic Dynamics & Control, 2010.

	Gary Anderson and George Moore. "An Efficient Procedure for Solving Linear Perfect Foresight Models."  Unpublished manuscript, Board of Governors of the Federal Reserve System. 1983.

}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
	\code{\link{callAMA}}
}


% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.


\examples{
hmat<-genHmat(system.file("extdata/example7.mod",package="AMA"),system.file("extdata/example7params.prm",package="AMA"))
bmat<-genBmat(hmat,4,1,1)
}


\keyword{ "models" }
\keyword{ "& algebra"}% __ONLY ONE__ keyword per line





