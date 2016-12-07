\name{candidateReduction}
\alias{candidateReduction}
\title{Ranks candidate variables for dimension reduction prior to
  estimation selection with the DSA} 
\description{
  Ranks explanatory variables based on their univariate association with
  an independent variable.  
}
\usage{
candidateReduction(formula, data, id = 1:nrow(data), family = gaussian,
weights=NULL, ...)
}
\arguments{
  \item{formula}{
    a symbolic description of the independent variable on which the
    dimension reduction is based: \code{formula} is typically set to 
    \code{Y ~ .} or \code{cbind(S,F) ~ .}. Currently supported outcomes are
    continuous or binomial (binary with 0s and 1s or a two-column
    matrix of successes in the first column and failures in the
    second).
  }
  \item{data}{
    a non-optional data frame containing both the response variable(s) as well as the
    candidate covariates. All variables in \code{data} will be
    considered candidate variables except for the outcome variable(s).
  }
  \item{id}{
    a vector identifying each independent experimental unit
    in the data with a unique value, i.e. repeated measurements are
    allowed in \code{data}. The length of \code{id} must correspond to the number of
    observations (\code{nrow(data)}). The default value for \code{id} is \code{1:nrow(data)} which
    indicates that all observations are independent.
  }
  \item{family}{
    currently either \code{binomial}, \code{multinomial} or \code{gaussian} (default). Used to determine whether
    logistic, multinomial (logit link only) or general linear models should be considered. 
  }
  \item{weights}{
    a vector of real numbers whose  number of elements is
    the number of observations. This vector contains
    the weights to be applied to each observation in \code{data}. The
    argument \code{weights} is ignored if the value for \code{weights} is
    \code{NULL} (default).
  }
  \item{...}{
    currently used internally to call \code{candidateReduction} within
    a \code{DSA} call. 
  }  
}
\details{
  The \code{candidateReduction} routine implements a dimension reduction
  step for the purpose of estimator selection. The candidate/explanatory
  variables are ranked based on their univariate association with an
  independent variable. Univariate regressions of the independent
  variable on each candidate variable are implemented and the p-values
  corresponding with the test of the null hypothesis (no association) are
  used to rank the candidate variables. Note that each model considered
  is fit based on semiparametric estimators. With repeated measurements,
  the working model for the matrix of variance-covariance for the
  residuals is diagonal
  (independence structure). The identity matrix is used if
  \code{family=gaussian} or \code{family=multinomial}. The inverse of
  the variance of the residuals is used in the diagonal if \code{family=binomial}. 

  This routine is called within the \code{DSA} routine for dimension
  reduction based on the association of each candidate variable with
  the outcome of interest. This routine can also be called externally by
  the user to further reduce the dimensionality of an estimator
  selection problem based on, for example, the association of each
  candidate variable with an exposure or treatment variable forced in
  the based model supplied to the \code{DSA}.
}
\value{
  a two-column matrix is returned. Each row contains the z statistics
  (for binomial or gaussian outcomes) or chi square
  statistics (for multinomial outcomes) and the corresponding p-value describing the association of the
  independent variable with a candidate variable identified by the row
  name. The matrix is ranked by increasing p-values. 
}
\author{
  Romain Neugebauer
}
\seealso{
  \code{\link{DSA}}
}
\examples{
library(DSA)

n <- 1000
W <- cbind(rnorm(n), rnorm(n) < 1, rnorm(n) < 2, rnorm(n, 2, 4), runif(n),
           rcauchy(n), rlogis(n) < .1, rnorm(n) < .1, rnorm(n, 120, 10),
           rnorm(n, 66, 2))

Y <- 10 + .5*W[,1] + .02*W[,1]^2 + .01*W[,1]*W[,2] + 2*W[,3] + .7*W[,4]^2 
Y <- as.matrix(as.integer(Y - mean(Y))/sd(Y))
trials <- rpois(n, lambda = 20)
successes <- sapply(1:n, function(i) {
  rbinom(1, size = trials[i], prob = pnorm(Y[i]))
})
failures <- trials - successes
colnames(W) <- paste("W", 1:ncol(W), sep = "")
data <- as.data.frame(cbind(W, "successes" = successes, "failures" = failures))

res <- candidateReduction(cbind(successes, failures) ~ ., data,family=binomial)  
res

res <- candidateReduction(W1 ~ ., as.data.frame(W),)  
res
}
\keyword{models}


















