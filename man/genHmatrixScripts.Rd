\name{genHmatrixScripts}
\alias{genHmatrixScripts}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Structural Coefficients Matrix Scripts}
\description{
	Creates an R file in the present working directory which, when run, creates the matrix H.
}
\usage{
genHmatrixScripts(modelFileNameFull)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{modelFileNameFull}{
String object containing the name of the model file.  Model file must be in MODELEZ syntax.
}

}
\details{
	Model input file must be in MODELEZ syntax.
}
\value{
	An R script, in the current directory, called "<model name>_sparseAimMatrices.r", which when run generates the matrix H.
}
\references{
For the details of MODELEZ syntax, see for example \url{http://www.federalreserve.gov/pubs/oss/oss4/papers/habitMatlab/habitMatlab.html}



	Gary S. Anderson. "A reliable and computationally efficient algorithm for imposing the saddle point property in dynamic models." Journal of Economic Dynamics & Control, 2010.

	Gary Anderson and George Moore. "An Efficient Procedure for Solving Linear Perfect Foresight Models."  Unpublished manuscript, Board of Governors of the Federal Reserve System. 1983.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
	\code{\link{callAMA}}
}

\examples{
genHmatrixScripts(system.file("extdata/example7.mod",package="AMA"))
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ "models" }
\keyword{ "& algebra"}% __ONLY ONE__ keyword per line




