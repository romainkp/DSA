\name{crossValidate}
\alias{crossValidate}
\title{Computation of the average cross-validated risk for any possible base model in the
  DSA routine}

\description{
  Calculate the average cross-validated risk based on the
  L2 loss-function for a base model that can be handled by the DSA routine .
}
\usage{
crossValidate(formula, data, ...)
}

\arguments{
  \item{formula}{a formula for a based model that can be handled by the DSA routine.}
  \item{data}{a non-optional data frame containing all variables in \code{formula}.}
  \item{\dots}{other arguments which are accepted by the DSA routine.}
}
\value{
A real number which corresponds to the average cross-validated risk associated
with the L2 loss function and the model specified with
\code{formula}. The data splits are defined by the DSA routine which is
used to compute the average cross-validated risk.
}
\seealso{\code{\link{DSA}}}
\examples{
library(DSA)

## an example using the state census data. (gaussian)  
data(state)
state.data <- as.data.frame(state.x77)
colnames(state.data) <- unlist(lapply(strsplit(colnames(state.data), " "),
                                      function(x) paste(x, collapse = "")))
res <- crossValidate(Murder ~ Population + Income + Illiteracy, data = state.data)
res

## an example with factors.
n <- 1000
inc <- round(runif(n,1,4))
edu <- round(runif(n,1,4))

W <- data.frame(income =
                as.factor(c("10000", "10000.35000", "35000.100000", "100000")[inc]),
                age = round(runif(n, 20, 60)),
                weight = rnorm(n, 140, 30),
                degree =
                as.factor(c("highschool", "highschool.grad", "college", "post.graduate")[edu]),
                male = rbinom(10, size = 1, prob = .5))

Vertical.Leap <- 8 + ifelse(inc == 1, 2, 0) + (60 - W$age) * .1 + W$male * 4 + rnorm(n, 0, 4)
data <- cbind(W, Vertical.Leap)                

res <- crossValidate(Vertical.Leap ~ income + degree, data = data, userseed=49,vfold=7,nsplits=2,silent=FALSE)
res
}
\keyword{models}


