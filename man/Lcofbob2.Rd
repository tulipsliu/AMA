\name{Lcofbob2.cof}
\alias{Lcofbob2.cof}
\docType{data}
\title{
Lcofbob2.cof Example Input data
}
\description{
Input matrix for example
}
\usage{Lcofbob2.cof}
\format{
Matrix elements columnwise.  }
\examples{
	hmat = read.table(system.file("extdata/Lcofbob2.cof",package="AMA"))
	hrows = 4
	hcols = 12
	neq = 4
	leads = 1
	lags = 1
	output <- callAMA(hmat, neq, leads, lags)
}
\keyword{datasets}
