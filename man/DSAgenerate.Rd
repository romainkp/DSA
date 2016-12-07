\name{DSAgenerate}
\alias{DSAgenerate}
\title{D/S/A algorithm to generate models without cross-validation}
\description{
  Generates models with the Deletion/Substitution/Addition (D/S/A) algorithm under user-specified
  constraints.
}
\usage{
DSAgenerate(formula, data, family = gaussian, weights = NULL, maxsize,
orderint, maxsumofpow, silent = TRUE, ...)
}
\arguments{
  \item{formula}{
    a symbolic description of the base model which specifies the
    independent/response variable(s) and all terms forced in
    the models returned. Typically, \code{formula} is set to
    \code{Y ~ 1} when no terms are forced in the models returned. Currently supported outcomes are
    continuous or binomial (binary with 0s and 1s or a two-column
    matrix of successes in the first column and failures in the
    second). Dependent/Explanatory candidate variables can be continuous
    or categorical (factors). Only polynomial main terms and
    polynomial interaction terms specified within the \code{I()}
    subroutine can be forced in the models returned, i.e. interactions
    between a factor and other variables are currently not suported.
  }
  \item{data}{
    a non-optional data frame containing both the response variable(s) as well as the
    candidate covariates to be considered in the model search.
  }
  \item{family}{
    currently either \code{binomial}, \code{multinomial} or \code{gaussian}. Used to determine whether
    logistic, multinomial (logit link only) or general linear models should be considered. 
  }
  \item{weights}{
    a vector of real numbers whose  number of elements is
    the number of observations. This vector contains
    the weights to be applied to each observation in the learning
    set (\code{data}). The argument \code{weights} is ignored if the value of
    \code{weights} is \code{NULL} (default).
  }
  \item{maxsize}{
    an integer strictly positive which limits the model search to
    candidate linear models with a maximum number of terms (excluding
    the intercept) lower or equal to
    \code{maxsize}. The argument \code{maxsize} must be larger or equal to the
    number of terms forced into the models through \code{formula}.
  }
  \item{orderint}{
    an integer strictly positive which limits the model search to
    candidate linear models with a maximum
    order of interactions equal to \code{orderint}. The argument
    \code{orderint} must be \code{1} if the number of terms forced in the models is
    equal to \code{maxsize}.
  }
  \item{maxsumofpow}{
    an integer larger or equal to \code{maxorderint} which limits the
    model search to candidate generalized linear models with terms involving
    variables whose powers sum up to a value lower or equal to
    \code{maxsumofpow}. The argument
    \code{maxsumofpow} must be \code{1} if the number of terms forced in the models is
    equal to \code{maxsize}.
  }
  \item{silent}{
    if \code{FALSE} then 
    intermediate messages will be printed to standard output showing the progress of the
    computations (message level set to 0). if \code{TRUE} then no
    intermediate messages will be printed to standard output (message
    level set to -1). One can further control the messaging level to be
    printed to standard output with a call to the \code{setDSAMessageLevel}
    subroutine preceding a call to the \code{DSA} routine where the argument
    \code{silent} is not referenced.
  }
  \item{...}{
    currently used internally to recursively call the DSA.
  }  
}
\details{
  The \code{DSAgenerate} routine implements the Deletion/Substitution/Addition
  (D/S/A) algorithm to generate candidate estimators defined with
  polynomial generalized linear models. Unlike the \code{DSA} routine,
  \code{DSAgenerate} does not perform data-adaptive etsimation based on
  cross-validation. Instead, \code{DSAgenerate} returns a list of
  estimators defined by models of size \code{1} to \code{maxsize} such
  that each estimator minimizes the empirical risk (not the
  cross-validated risk) among all estimators of the same size considered.
  The space of considered candidate estimators is parameterized with three variables:
  \code{maxsize}, \code{orderint} and \code{maxsumofpow}. 
  
  The D/S/A algorithm is an aggressive model search algorithm which
  iteratively generates polynomial generalized linear models based on
  the existing terms in the current 'best' model and the following three steps: 1) a deletion step which removes a term from the model,
  2) a substitution step which replaces one term with another, and 3) an
  addition step which adds a term to the model. The search for the 'best'
  estimators starts with the base model specified with \code{formula}:
  typically the intercept model except
  when the user requires a number of terms to be forced in the models returned.

  The search for the 'best' estimators is limited by three user-specified
  arguments: \code{maxsize}, \code{orderint} and \code{maxsumofpow}. The first argument
  limits the maximum number of terms in the models considered (excluding
  the intercept). The second
  argument limits the maximum order of interactions for the models
  considered. All terms in the models considered are composed of interactions of
  variables raised to a given power. The third argument limits the
  maximum sum of powers in each term.
}
\value{
  a list of four objects is returned: 

  \item{models.allsizes}{list of the 'best' models for each size
    considered and corresponding with the order of interaction equal to
    \code{orderint}.}
  \item{family}{a description of the link function defining the
    candidate models considered: \code{gaussian} indicates
    the indentity function while \code{binomial} and \code{multinomial} indicates  the logit function.}
  \item{computing.time}{an object of class \code{difftime} which contains the
    information about the computing time associated with this call to \code{DSAgenerate}.}
  \item{call}{the DSAgenerate call which generated this object.}
}
\author{
  Romain Neugebauer and James Bullard based on the original C code from
  Sandra Sinisi.
}
\seealso{
  \code{\link{DSA}}, \code{\link{setDSAMessageLevel}},
  \code{\link{getDSAMessageLevel}}, \code{\link{difftime}}.
}
\references{
  1. Sandra E. Sinisi and Mark J. van der Laan, "Loss-Based
  Cross-Validated Deletion/Substitution/Addition Algorithms in Estimation"
  (March 2004). U.C. Berkeley Division of Biostatistics Working Paper
  Series. Working Paper 143.
  http://www.bepress.com/ucbbiostat/paper143
}
\examples{
library(DSA)

##
## an example using binomial - two column outcome and forced terms
##
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
colnames(W) <- paste("V", 1:ncol(W), sep = "")
data <- as.data.frame(cbind(W, "successes" = successes, "failures" = failures))

res <- DSA(cbind(successes, failures) ~ 1, data = data, family = binomial, maxsize = 8,
           maxorderint = 2, maxsumofpow = 3)  
summary(res)
plot(res)

## one may prefer to use a model with an order of interaction equal to 2
res <- DSAgenerate(cbind(successes, failures) ~ 1, data = data, family = binomial, maxsize = 8,
           orderint = 2, maxsumofpow = 3)  
res

## now add some forced terms
res <- DSAgenerate(cbind(successes, failures) ~ V1 + I(V4^2), data = data, family = binomial, maxsize = 8,
           orderint = 2, maxsumofpow = 3) 
res
}

\keyword{models}


















