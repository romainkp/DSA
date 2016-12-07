\name{plot.DSA}
\alias{plot.DSA}
\title{A visual diagnostic tool for data-adaptive estimation with the
  DSA routine}
\description{
  Plot the average cross-validated risks obtained from a call to the DSA routine. 
}
\usage{
plot.DSA(x, plot.compare = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{an object of class 'DSA'. }
  \item{plot.compare}{if \code{TRUE}, plots are overlaid to
    facilitate the comparison of average cross-validated risks for
    all sizes, interaction orders and dimension reduction levels considered.}
  \item{\dots}{currently additional arguments are ignored. }
}
\value{
  A list of plots. Currently model sizes are plotted against average
  cross-validated risks for each dimension reduction level considered separately for each
  interaction order if \code{plot.compare=FALSE} or on a
  single plot if \code{plot.compare=TRUE}.
}
\author{James Bullard.}
\seealso{\code{\link{DSA}}.}
\examples{
library(DSA)

## an example using the state census data. (gaussian)  
data(state)
state.data <- as.data.frame(state.x77)
colnames(state.data) <- unlist(lapply(strsplit(colnames(state.data), " "),
                                      function(x) paste(x, collapse = "")))
res <- DSA(Murder ~ 1, data = state.data, maxsize = 5, maxsumofpow = 2, maxorderint = 2)
plot(res)
}
\keyword{models}% at least one, from doc/KEYWORDS

