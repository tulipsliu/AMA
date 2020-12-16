\name{example7params}
\alias{example7params}
\docType{data}
\title{
example7params.txt -- Text input for example
}
\description{
Contains parameter assignments for the example7 matrix assignments.
}
\format{
Rows of assignment statements for model parameters.
}


\examples{
hmat<-genHmat(system.file("extdata/example7.mod",package="AMA"),system.file("extdata/example7params.prm",package="AMA"))
bmat<-genBmat(hmat,4,1,1)
}
\keyword{datasets}
