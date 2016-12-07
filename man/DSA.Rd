\name{DSA}
\alias{DSA}
\title{Data-Adaptive Estimation with Cross-Validation and the D/S/A Algorithm}
\description{
  Performs data-adaptive estimation through estimator selection based on
  cross-validation and the L2
  loss function. Candidate estimators are defined with polynomial generalized linear models generated with
  the Deletion/Substitution/Addition (D/S/A) algorithm under user-specified
  constraints.
}
\usage{
DSA(formula, data, id = 1:nrow(data), family = gaussian, weights = NULL,
candidate.rank = NULL, rank.cutoffs =  NULL, maxsize, maxorderint,
maxsumofpow, Dmove=TRUE, Smove=TRUE, usersplits = NULL, userseed = NULL,
vfold = 5, nsplits = 1, model.matrices = FALSE, silent = TRUE, ...)
}
\arguments{
  \item{formula}{
    a symbolic description of the base model which specifies the
    independent/response variable(s) and all terms forced in
    the final model. Typically, \code{formula} is set to
    \code{Y ~ 1} when no terms are forced in the final model. Currently supported outcomes are
    continuous, binomial (binary with 0s and 1s or a two-column
    matrix of successes in the first column and failures in the
    second) and multinomial (categorical specified as factors). Dependent/Explanatory candidate variables can be continuous
    or categorical (factors). Only polynomial main terms and
    polynomial interaction terms specified within the \code{I()}
    subroutine can be forced in the final model, i.e. interactions
    between a factor and other variables are currently not suported.
  }
  \item{data}{
    a non-optional data frame containing both the response variable(s) as well as the
    candidate covariates to be considered in the data-adaptive
    estimation procedure.
  }
  \item{id}{
    a vector identifying each independent experimental unit
    in the data with a unique value, i.e. repeated measurements are
    allowed in \code{data}. The length of \code{id} must correspond to the number of
    observations (\code{nrow(data)}). The default value for \code{id} is \code{1:nrow(data)} which
    indicates that all observations are independent.
  }
  \item{family}{
    currently either \code{binomial}, \code{multinomial}, or \code{gaussian}. Used to determine whether
    logistic, multinomial logistic (softmax function), or general linear models should be considered. 
  }
  \item{weights}{
    a matrix of real numbers whose number of rows
    is the number of data splits plus one (i.e., \code{vfold+1} or \code{nrow(usersplits)+1}) and whose number of columns is
    the number of observations. The first \code{vfold} or \code{nrow(usersplits)} rows contain the weights to be applied to
    each observation (whether in the training or validation set) at each
    split of the data. The last row contains
    the weights to be applied to each observation in the learning
    set. The argument \code{weights} is ignored if the value of
    \code{weights} is \code{NULL} (default).
  }
  \item{candidate.rank}{a vector of named real numbers which ranks all
    candidate variables in \code{data}. Missing values are not allowed. Each element's name must
    correspond to the name of one of the candidate variables in
    \code{data}. The candidate variables' ranks are used to
    determine which set(s) of covariates in \code{data} will be
    considered in the DSA call based on the argument \code{rank.cutoffs}.  When
    \code{candidate.rank} is \code{NULL} (default), the variables in
    data are ranked with a call to \code{candidateReduction}. Note that
    the candidate dummy variables automatically created for all factors in \code{data} receive different
    ranks by default while their ranks are identical when \code{candidate.rank} is user-supplied. 
  }
  \item{rank.cutoffs}{a vector of real numbers which is used with
    \code{candidate.rank} to define the set(s) of candidate variables in
    \code{data} that will be considered in the data-adaptive estimation procedure. By default,
    \code{rank.cutoffs} is set to \code{NULL}, which means that all candidate
    variables supplied in \code{data} are considered in the DSA call.
  }
  \item{maxsize}{
    an integer strictly positive which limits the data-adaptive
    estimation procedure to candidate linear models with a maximum
    number of terms (excluding the intercept) lower or equal to
    \code{maxsize}. The argument \code{maxsize} must be larger or equal to the
    number of terms forced into the model through \code{formula}.
  }
  \item{maxorderint}{
    an integer strictly positive which limits the data-adaptive
    estimation procedure to candidate linear models with a maximum
    order of interactions lower or equal to \code{maxorderint}.
  }
  \item{maxsumofpow}{
    an integer larger or equal to \code{maxorderint} which limits the data-adaptive
    estimation procedure to candidate generalized linear models with terms involving
    variables whose powers sum up to a value lower or equal to
    \code{maxsumofpow}.
  }
  \item{Dmove}{
    a logical parameter which specifies whether deletion moves are used to search through
    the space of candidate estimators. By default, deletion moves are allowed.
  }
  \item{Smove}{
    a logical parameter which specifies whether substitution moves are used to search through
    the space of candidate estimators. By default, substitution moves are allowed.
  }  
   \item{usersplits}{
    an optional matrix of 0s or 1s whose number of columns is the number of
    observations. Each row indicates how the learning set
    is split in a training (0s) and validation (1s) set. The default
    value for \code{usersplits} is \code{NULL} which indicates that the data should be
    randomely split based on the v-fold splitting scheme corresponding to the values
    for \code{userseed}, \code{vfold}, \code{nsplits}, and \code{id}.
  }
  \item{userseed}{
    a single value interpreted as an integer and used to
    set the seed of the random number
    generator which determines the v-fold data splits.  If \code{userseed=NULL} (default) then
    the current seed from the current R session is used. The argument
    \code{userseed} is ignored if \code{usersplits} is non-null.}
  \item{vfold}{
    an integer strictly larger than 1 specifying the number of splits
    (i.e, v) in
    the default v-fold splitting scheme for the
    cross-validation procedure. The default value for \code{vfold} is
    \code{5}. The argument \code{vfold} is ignored if the argument
    \code{usersplits} is non-null.
  }
  \item{nsplits}{
    an integer larger than 1 indicating the number of v-fold splits for the
    cross-validation procedure. The default value for \code{nsplits} is
    \code{1}. This argument is ignored if the user provides the
    data splitting scheme with \code{usersplits}.
  }
  \item{model.matrices}{
    a logical variable indicating that only the design and response
    matrices should be returned.
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
  The DSA routine implements a general data-adaptive estimation
  procedure based on cross-validation and the L2 loss function. The final estimator is selected from a set
  of candidate estimators defined with polynomial generalized linear models generated by the
  Deletion/Substitution/Addition (D/S/A) algorithm. The space of
  candidate estimators is parameterized with three or four variables:
  \code{maxsize}, \code{maxorderint}, \code{maxsumofpow} and
  \code{rank.cutoffs} (optional). The final
  model returned minimizes the empirical risk on the
  learning set among all estimators considered and characterized by the
  "optimum" size, order of interactions and set of candidate variables selected by cross-validation.
  
  The D/S/A algorithm is an aggressive model search algorithm which
  iteratively generates polynomial generalized linear models based on
  the existing terms in the current 'best' model and the following three steps: 1) a deletion step which removes a term from the model,
  2) a substitution step which replaces one term with another, and 3) an
  addition step which adds a term to the model. Note that the user can inhibit deletion and substitution moves independenlty with \code{Dmove} and \code{Smove}. The search for the 'best'
  estimator starts with the base model specified with \code{formula}:
  typically the intercept model except
  when the user requires a number of terms to be forced in the final model.

  The search for the 'best' estimator is limited by three or four user-specified
  arguments: \code{maxsize}, \code{maxorderint}, \code{maxsumofpow} and \code{rank.cutoffs} (optional). The first argument
  limits the maximum number of terms in the models considered (excluding
  the intercept). The second
  argument limits the maximum order of interactions for the models
  considered. All terms in the models considered are composed of interactions of
  variables raised to a given power. The third argument limits the
  maximum sum of powers in each term. The fourth argument limits the set of candidate variables to
  be considered in each model. Only the variables whose ranks specified
  by \code{candidate.rank} are smaller than the threshold(s) specified by
  \code{rank.cutoffs} are allowed in the models considered (and, if applicable, the variable(s) forced in the model). Note that
  the default ranking of the candidate variables in \code{data} is based
  on their univariate association with the independent variable(s) and
  obtained with the routine \code{candidateReduction} if \code{candidate.rank} is \code{NULL}. If \code{rank.cutoffs} is \code{NULL} then all candidate variables in \code{data} are allowed in the models considered.

  This data-adaptive estimation procedure allows comparison of models
  based on different number of observations and can account for
  informative censoring through the use of weights in each
  regression. These weights are provided to the DSA routine with the
  argument \code{weights}.

  The DSA routine currently supports data-adaptive estimation
  for continuous, binomial and multinomial outcomes. When the outcome is binomial, the estimators considered are
  based on polynomial generalized linear models where the link function is the
  logit function. When the outcome is categorical, the estimators considered are
  based on multinomial logit models. Factors can be candidate variables with the caveat
  that there are currently limitations
  to the use of factors in terms forced in the final model (see \code{formula} above). 

  The default cross-validation splitting scheme  is v-fold
  where the value for v is specified with the argument \code{vfold}. The DSA
  routine performs the data splits based on the value for
  \code{vfold}, \code{nsplits}, \code{id} and the \code{userseed}
  arguments. The argument \code{nsplits} specifies the number of v-fold
  splits, e.g. if \code{nsplits=2} and \code{vfold=5} then the data is split twice based on
  the 5-fold splitting scheme and thus the DSA call relies on 10 data splits. The
  argument \code{id}
  identifies the independent experimental units in the data and ensures
  that the training and validation sets are independent. 
  The argument \code{userseed} is
  used to set the seed of the R random number generators with the routine
  \code{set.seed()}. This allows for reproducible results. The user can specify
  an alternative cross-validation splitting scheme with the argument \code{usersplits}. 
}
\value{
  If \code{model.matrices=TRUE} then a list of two objects is returned, \code{X} and \code{Y}, corresponding respectively with the
  matrices of all candidate variables and outcome in \code{data}.
  Otherwise, an object of class 'DSA' is returned. The traditional \code{print},
  \code{summary}, \code{predict}, and \code{plot} methods are then available on the returned
  object. Objects of class 'DSA' have the following attributes: 
  
  \item{average.CVrisks}{a matrix or tridimensional array containing the average cross-validated risks
    defined by the L2 loss function over the splits of the data for estimators (models) indexed by
    1) a size ranging from 1 to \code{maxsize+1} (intercept included);
    2) a maximum order of interactions ranging from 1 to \code{maxorderint}; and
    3) by the dimension reduction levels specified in \code{rank.cutoffs} (if non-null).
    The average cross-validated risk with coordinates (i,j,k) corresponds to the average cross-validated risk
    for estimators (models) indexed by  a maximum order of interactions
    i, size j (intercept included), and a rank cut-off value identified
    by the column name.
    Note that the average cross-validated risk for a given size, order of interactions and dimension reduction level is the average across data splits of the
    average risks computed on the validation sets. The average
    risks on a given validation set corresponds to the residuals sum of
    squares (missing residuals excluded) divided by the
    number of (complete) observations in the validation set.}
  \item{min.risk}{the coordinates identifying the minimum average cross-validated risk in
    \code{average.CVrisks}. Note that the corresponding average cross-validated
    risk does not typically correspond to the average cross-validated risk for the final
  model selected by the DSA routine.}
   \item{final.cutoff}{the 'best' or user-specified dimension reduction
     level corresponding with the model selected.}
\item{model.selected}{a formula representing the final model selected
    by the DSA routine. The variables in \code{model.selected}
    correspond to variables
    in \code{data}.}
  \item{family}{a description of the link function for the final
    generalized linear model selected: \code{gaussian} indicates
    the indentity function while \code{binomial} and \code{multinomial} indicates  the logit function.}
  \item{coefficients}{the coefficients of the final model fitted on the
    learning set.}
  \item{average.CVrisk.model.selected}{the average cross-validated risk
    for \code{model.selected}.}  
  \item{models.allsizes}{list of the 'best' models of each size considered on the learning
  set and corresponding with the 'optimum' order of interaction and dimension reduction level selected by
  cross-validation.}
  \item{candidate.vars}{list of the candidate variables considered in
    the data-adaptive estimation procedure.
    Note that this list does not necessarily contain all candidate variables in \code{data} (except when \code{rank.cutoffs} is \code{NULL}) but only contain the variables whose ranks are lower than or equal to the largest cut-off
    value in \code{rank.cutoffs}. Note that this list includes the variables forced in the final
    model, if applicable.}
  \item{selected.vars}{list of the variables corresponding with the
    'best' dimension reduction level, i.e. all variables whose ranks are
    lower than or equal to the 'best' cut-off value selected by
    cross-validation. Note that this list includes the variables forced in the final
    model, if applicable. If \code{rank.cutoffs} is \code{NULL} or
    contains only one valid cut-off value then selected.vars is \code{NULL}.}
  \item{candidate.rank}{the original argument \code{candidate.rank}
    passed in to the DSA routine or, if not specified, the default ranking of all variables in \code{data} obtained with \code{candidateReduction}. If \code{rank.cutoffs} is \code{NULL} then \code{candidate.rank} should be ignored.
  }
  \item{rank.cutoffs}{the original argument \code{rank.cutoff}
    passed in to the DSA routine or, if not-specified, its default value \code{NULL}.
    }
  \item{computing.time}{an object of class \code{difftime} which contains the
    information about the computing time associated with this DSA object.}
  \item{call}{the DSA call which generated this DSA object.}
  \item{formula}{the original base formula passed in to the DSA routine.}
  \item{data}{the original data frame passed in to the DSA routine. }
  \item{id}{the value of \code{id} used in the DSA call.}
  \item{userseed}{\code{NULL} if the argument \code{usersplits} was provided by the user or
    if the user did not provide \code{userseed}, otherwise equal to \code{userseed}. }
  \item{splits}{a matrix or a list of matrices that identifies how
    the data were split, i.e. the \code{usersplits} matrix if provided by the user or
    the corresponding matrix (when \code{nsplits = 1}) or list of matrices
    (when \code{nsplits > 1}) randomely generated by the DSA routine.}
  \item{weights}{the weights used in the original call.}
  \item{moves}{a logical vector indicating whether deletion and
    substitution moves where permitted or inhibited in the model search.}
}
\note{
  The computing time for a DSA call can be very large depending on the input arguments which define the intensity of the estimator selection. DSA calls are increasingly more computing intensive with binomial and multinomial outcomes. It is often recommended to run this routine in BATCH mode.
}
\author{
  Romain Neugebauer and James Bullard based on the original C code from
  Sandra Sinisi.
}
\seealso{
  \code{\link{I}}, \code{\link{set.seed}},
  \code{\link{difftime}}, \code{\link{setDSAMessageLevel}},
  \code{\link{getDSAMessageLevel}}, \code{\link{crossValidate}},
  \code{\link{candidateReduction}}, \code{\link{DSAgenerate}},
  \code{\link{summary.DSA}},\code{\link{coefficients.DSA}},
  \code{\link{residuals.DSA}}, \code{\link{predict.DSA}},
  \code{\link{plot.DSA}}.
}
\references{
  1. Mark J. van der Laan and Sandrine Dudoit, "Unified Cross-Validation
  Methodology For Selection Among Estimators and a General Cross-Validated
  Adaptive Epsilon-Net Estimator: Finite Sample Oracle Inequalities and
  Examples" (November 2003). U.C. Berkeley Division of Biostatistics
  Working Paper Series. Working Paper 130.
  http://www.bepress.com/ucbbiostat/paper130

  2. Mark J. van der Laan, Sandrine Dudoit, and Aad W. van der Vaart, "The
  Cross-Validated Adaptive Epsilon-Net Estimator" (February 2004). U.C.
  Berkeley Division of Biostatistics Working Paper Series. Working Paper
  142.
  http://www.bepress.com/ucbbiostat/paper142
  
  3. Sandra E. Sinisi and Mark J. van der Laan, "Loss-Based
  Cross-Validated Deletion/Substitution/Addition Algorithms in Estimation"
  (March 2004). U.C. Berkeley Division of Biostatistics Working Paper
  Series. Working Paper 143.
  http://www.bepress.com/ucbbiostat/paper143
}
\examples{
library(DSA)

## an example using the state census data. (gaussian)  
data(state)
state.data <- as.data.frame(state.x77)
colnames(state.data) <- unlist(lapply(strsplit(colnames(state.data), " "),
                                      function(x) paste(x, collapse = "")))
res <- DSA(Murder ~ 1, data = state.data, maxsize = 5, maxsumofpow = 2, maxorderint = 2)
residuals(res)
coefficients(res)
summary(res)

res <- DSA(Murder ~ 1, data = state.data, maxsize = 5, maxsumofpow = 2, maxorderint = 2,nsplits=10,silent=FALSE)
residuals(res)
coefficients(res)
summary(res)

## an example using weights (gaussian).
ddir <- paste(system.file(package = "DSA"), "/", "data", "/", sep = "")
x <- read.table(paste(ddir, "smalllymphRoX.txt", sep = ""), nrow=240)
y <- read.table(paste(ddir, "logT1lymph.txt", sep = ""))
weights <- as.matrix(read.table(paste(ddir, "lymphWeights_111505", sep = "")))
colnames(y) <- "loglymph"
data <- cbind(x, y)

res <- DSA(loglymph ~ 1, data = data, maxsize = 5, maxorderint = 2, maxsumofpow = 2,
           weights = weights)
plot(res)
summary(res)

## an example using simulated data (binomial, logistic link).
n <- 1000
W <- cbind(rnorm(n), rnorm(n) < 1, rnorm(n) < 2, rnorm(n, 2, 4), runif(n),
           rcauchy(n), rlogis(n) < .1, rnorm(n) < .1, rnorm(n, 120, 10),
           rnorm(n, 66, 2))

Y <- 22 + .5*W[,1] + .02*W[,1]^2 + .01*W[,1]*W[,2] + 2*W[,3] + .7*W[,4]^2 
Y <- as.matrix(as.integer(Y - mean(Y))/sd(Y))
Y <- as.matrix(as.integer(Y < rlogis(n))) 
colnames(W) <- paste("V", 1:ncol(W), sep = "")
colnames(Y) <- "Y"
data <- as.data.frame(cbind(W, Y))

res <- DSA(Y ~ 1, data = data, family = binomial, maxsize = 8, maxorderint = 2, maxsumofpow = 3)  
plot(res)
summary(res)


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
res <- DSA(Vertical.Leap ~ 1, data = data, maxsize = 8, maxorderint = 2, maxsumofpow = 2)  
plot(res, plot.compare = TRUE)
summary(res)

## now add some forced terms.
res <- DSA(Vertical.Leap ~ degree, data = data, maxsize = 8, maxorderint = 2, maxsumofpow = 2)  
plot(res, plot.compare = TRUE)
summary(res)


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


## now add some forced terms
res <- DSA(cbind(successes, failures) ~ V1 + I(V4^2), data = data, family = binomial, maxsize = 8,
           maxorderint = 2, maxsumofpow = 3) 
summary(res)
plot(res)

## with user-specified dimension reduction
res <- DSA(cbind(successes, failures) ~ V1, data=data,family=binomial,  rank.cutoffs=0.1, maxsize = 5, maxorderint = 1, maxsumofpow = 2,userseed=29,nsplits=2,silent=FALSE)
summary(res)
res$candidate.vars
res$selected.vars

## with dimension reduction based on cross-validation
res <- DSA(cbind(successes, failures) ~ V1, data=data,family=binomial,  rank.cutoffs=c(0.1,0.2,0.6), maxsize = 5, maxorderint = 1, maxsumofpow = 2,userseed=29)
summary(res)
res$candidate.vars
res$selected.vars

## now a multinomial model.
NN <- 1000
W <- cbind(rep(1, NN), runif(NN, 0, 5), rnorm(NN, 5, 2))

PY <- apply(W, 1, function(x) {
coefs <- matrix(c(0,0,0,-9, 2, .8,-6, 1, .7,-4, .5, .6), nrow = 4, ncol = 3, byrow = TRUE)
denom <- sum(apply(coefs, 1, function(cc) {exp(x \%*\% cc)}))
as.numeric(apply(coefs, 1, function(cc) {exp(x \%*\% cc)/denom}))
})

## add some columns to W which don't matter
W <- cbind(W, runif(NN, 0, 5), rnorm(NN, 5, 2))

## turn them into counts.
Y <- apply(PY, 2, function(x) {
  ru <- runif(1)
  if (ru <= x[1])
    0
  else if (ru <= x[2] + x[3])
    1
  else if (ru <= x[3] + x[4])
    2
  else
    3
})

dta.2 <- data.frame(W[,-1], "Y" = factor(Y))
levels(dta.2[,"Y"]) <- c("B","A","D","C")
res <- DSA(Y ~ 1, family=multinomial, data = dta.2, maxsize = 5,
           maxorderint = 2, maxsumofpow = 2,vfold = 5, userseed = 87,
           model.matrices = FALSE, nsplits=1)
res
summary(res)

## Another example with the state census data where deletion and/or substitution moves are inhibited
res1 <- DSA(Murder ~ 1, data = state.data, maxsize = 10, maxsumofpow = 2, maxorderint = 2)
summary(res1)

res2 <- DSA(Murder ~ 1, data = state.data, maxsize = 10, maxsumofpow = 2, maxorderint = 2,Dmove=FALSE)
summary(res2)

res3 <- DSA(Murder ~ 1, data = state.data, maxsize = 10, maxsumofpow = 2, maxorderint = 2,Smove=FALSE)
summary(res3)

res4 <- DSA(Murder ~ 1, data = state.data, maxsize = 10, maxsumofpow = 2, maxorderint = 2,Dmove=FALSE,Smove=FALSE)
summary(res4)
}

\keyword{models}


















