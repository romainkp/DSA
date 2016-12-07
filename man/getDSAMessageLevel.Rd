\name{getDSAMessageLevel}
\alias{getDSAMessageLevel}
\title{Get current messaging level for the DSA routine}
\description{
Get the current complexity level for intermediate messages to be printed to
standard output during a call to the DSA routine.
}
\usage{
getDSAMessageLevel()
}
\value{
The current value (integer) for the messaging level.
}
\author{James Bullard.}
\seealso{\code{\link{DSA}}.}
\examples{
library(DSA)

## an example using the state census data. (gaussian)  
data(state)
state.data <- as.data.frame(state.x77)
colnames(state.data) <- unlist(lapply(strsplit(colnames(state.data), " "),
                                      function(x) paste(x, collapse ="")))

getDSAMessageLevel()
res <- DSA(Murder ~ 1, data = state.data, maxsize = 5, maxsumofpow = 2,
           maxorderint = 2, userseed=78)
res

setDSAMessageLevel(0) #default value
getDSAMessageLevel()
res <- DSA(Murder ~ 1, data = state.data, maxsize = 5, maxsumofpow = 2,
           maxorderint = 2 ,userseed=78)
res
}
\keyword{models}
