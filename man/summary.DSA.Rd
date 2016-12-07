\name{summary.DSA}
\alias{summary.DSA}
\title{Summarization of the results of a call to the DSA routine}
\description{
  Provides a useful summary of a call to the DSA routine.
}
\usage{
summary.DSA(object, ...)
}


\arguments{
  \item{object}{an object of class 'DSA'.}
  \item{\dots}{currently ignored.}
}

\value{An object of class 'summary.DSA'. Objects of class 'summary.DSA'
  have the following attributes:
  \item{model.selected}{a formula representing the final model selected
    by the DSA routine.}
  \item{family}{a description of the link function for \code{model.selected}: \code{gaussian} indicates
    the indentity function while \code{binomial} indicates  the logit function.}
  \item{n.obs}{number of observations in the data frame \code{data} passed in to the DSA routine.}
  \item{n.nonna}{number of complete observations on which \code{model.selected}
    is based (no missing values).}
  \item{best.average.CVrisk}{minimum average cross-validated risk among
    all average cross-validated risks
    associated with all model sizes, interaction orders and dimension
    reduction levels considered by
  the DSA routine.}
  \item{coefficients}{the coefficients of \code{model.selected} fitted on
    \code{n.nonna} observations from the learning set.}
  \item{io.and.size}{the 'best' model size and order of interaction
    selected by cross-validation}  
  \item{average.CVrisk.model.selected}{average cross-validated risk
    for \code{model.selected}.}
  \item{n.id}{number of independent experimental units in the data frame
    \code{data} passed in to the DSA routine.}
  \item{n.id.nonna}{number of complete independent experimental units on which \code{model.selected}
    is based (no missing values).}
  \item{nselected.vars}{the number of candidate variables after
    dimension reduction based on a user-specified rank cut-off or a rank
  cut-off selected by cross-validation.}
  \item{ncandidate.vars}{the number of candidate variables truly
    considered, i.e. the number of variables whose ranks are lower than or equal to the largest cut-off
    value in \code{rank.cutoffs}.}
  \item{n.vars}{the number of candidate variables considered
    pre-dimension reduction.}
  \item{rank.cutoffs}{the rank cut-offs considered.}
  \item{best.rank.cutoff}{\code{NULL} if the dimension reduction level
    is user-specified or the 'best' dimension reduction level if the
    dimension reduction step is based on cross-validation.}
  \item{formula}{the original base formula passed in to the DSA routine.}
  \item{moves}{a logical vector indicating whether deletion and
    substitution moves where permitted or inhibited in the model search.}
}
\author{
  James Bullard.
}
\seealso{
  \code{\link{DSA}}.
}
\examples{
library(DSA)

## an example using the state census data. (gaussian)  
data(state)
state.data <- as.data.frame(state.x77)
colnames(state.data) <- unlist(lapply(strsplit(colnames(state.data), " "),
                                      function(x) paste(x, collapse = "")))
res <- DSA(Murder ~ 1, data = state.data, maxsize = 5, maxsumofpow = 2,maxorderint = 2)
res
summary(res)
}

\keyword{models}

