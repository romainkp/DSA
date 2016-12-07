\name{setDSAMessageLevel}
\alias{setDSAMessageLevel}
\title{Control of the messaging level for the DSA routine}
\description{
Set the complexity level for intermediate messages to be printed to
standard output during a call to the DSA routine.
}
\usage{
setDSAMessageLevel(newlevel)
}
\arguments{
  \item{newlevel}{an integer ranging from -1 to 4. the complexity
    of the messaging level increases as the value for \code{newlevel}
    gets larger.}
}
\value{
The passed in value for \code{newlevel}.
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

setDSAMessageLevel(0)
res <- DSA(Murder ~ 1, data = state.data, maxsize = 5, maxsumofpow = 2,
           maxorderint = 2, userseed=78)
res

setDSAMessageLevel(-1) #default value
res <- DSA(Murder ~ 1, data = state.data, maxsize = 5, maxsumofpow = 2,
           maxorderint = 2 ,userseed=78)
res
}
\keyword{models}

