##
## public method (s)
##

set.DSA.fit.debug <- function(debug.level = 0) {
  invisible(.C("R_DSA_set_debug_level", as.integer(debug.level)))
}

##
## non-public methods
##

.DSA.transpose <- function(x) {
  storage.mode(x) <- "double"
  z = matrix(double(nrow(x)*ncol(x)), nrow = nrow(x))
  
  .C("R_DSA_transpose", x, nrow(x), ncol(x), z)
}


##
## performs logistic regression given a
## design matrix and a vector of outcomes. 
##
.DSA.do.lr <- function(x, y, p = NULL, weights = rep(1, nrow(x)),
                       max.iter = 25, transpose = TRUE) {
  n   <- as.integer(nrow(x))
  tdx <- as.integer(ncol(x))
  ny  <- as.integer(1)
  tol <- as.double(1e-8)

  if (is.null(p))
    p <- as.integer(tdx)
  else
    p <- as.integer(p)

  coefficients <- double(p)
  residuals <- double(n)
  effects <- double(n)
  rank <- integer(1)
  pivot <- 1:p 
  qraux <- double(p)
  dev <- double(1)
  
  if (is.null(weights)) {
    weights <- as.double(rep(1, n))
  }
  else {
    weights <- as.double(weights)
  }

  max.iter <- as.integer(max.iter)
  new.coeffs <- double(p)
  xtwx <- matrix(double(n*tdx), nrow = n)
  dim(xtwx) <- dim(x)
  yshift <- double(n)
  work <- double(2 * p)

  # column based storage. 
  x <- t(x)
  
  .C("R_DSA_fit_glm_lr", qr = as.double(x), n, tdx, p, y = as.double(y), tol, coefficients = coefficients,
     residuals, effects, rank, pivot, qraux, work = work,
     dev = dev, weights = weights, max.iter, new.coefficients = new.coeffs, xw = xtwx,
     yw = yshift, do.transpose = as.integer(0), fail = integer(1))
}
                  

##
## performs ordinary least squares.
##
.DSA.do.ols <- function(x, y, p = NULL, weights = NULL, my.tol = 1e-6) {
  n   <- as.integer(nrow(x))
  tdx <- as.integer(ncol(x))
  tol <- as.double(my.tol)

  if (is.null(p)) {
    p <- as.integer(tdx)
  }
  else {
    p <- as.integer(p)
  }
  
  if (is.null(weights)) {
    weights <- as.double(rep(1, n))
  }
  else {
    weights <- as.double(weights)
  }
  
  coefficients <- double(p)
  residuals <- double(n)

  effects <- double(n)
  rank <- integer(1)
  pivot <- as.integer(1:p) 
  qraux <- double(p)

  # these routines use column based storage. 
  x <- t(x)
 
  .C("R_DSA_fit_mlr", qr = as.double(x), n = n, tdx = tdx, ip = p, y = as.double(y),
     weights = weights, tol = as.double(tol), coefficients = coefficients, 
     residuals = residuals, effects = effects, rank = integer(1), pivot = pivot, 
     qraux = qraux, work = double(2 * p), rss = double(1),
     do.transpose = as.integer(0), fail = integer(1))
}

##
## generate a design matrix to test the model fitting
## routines. 
##
.DSA.generate.X <- function(n, outcome = 'binary') {
  x <- cbind(rep(1,n), (115 + ((1:n)/(n/20)) + rnorm(n)),
             round(6 + (1:n)/(n/10) + rnorm(n, 0, 2)), (60 + (1:n)/n + rnorm(n, 0, 6)),
             rnorm(n) > .5)

  x <- cbind(x, x[,2]^2 + rnorm(n), x[,2]*x[,3] + rnorm(n), x[,2]*x[,4] + rnorm(n),
             x[,3]*x[,4] + rnorm(n), x[,3]^2 + rnorm(n))
 
  
  
  rs <- c(0, .0001, .0001, .0001, .6, .2, .01, .04, .05, .002)
  
  theta <- x %*% rs
  
  if (outcome == 'continuous')
    y <- theta + rnorm(n, 0, 1)
  else
    y <- rlogis(n) < (theta - mean(theta))/sd(theta)

  
  
  list(x = x, y = y, theta = theta)
} 


##
## generate candidate matrix for the DSA to perform model selection. 
##
.DSA.generate.W <- function(n = 100, outcome = 'continuous') {
  W <- cbind(rnorm(n), rnorm(n) < 1, rnorm(n) < 2, rnorm(n, 2, 4), runif(n),
             rcauchy(n), rlogis(n) < .1, rnorm(n) < .1, rnorm(n, 120, 10),
             rnorm(n, 66, 2))

  Y <- as.matrix((22 + .5*W[,1] + .02*W[,1]^2 + .01*W[,1]*W[,2] + 2*W[,3] + .7*W[,4]^2 +
                  -1.2*W[,5] + 1.1*W[,6]*W[,7] + -.5*W[,8]^2 + .01*W[,9] + .001*W[,10]*W[,9]^2))

  if (outcome == 'binary') {
    Y <- rlogis(n) > (Y - mean(Y))/sd(Y) 
  }

  colnames(W) <- paste("V", 1:ncol(W), sep = "")
  colnames(Y) <- "Y"
  
  list(w = W, y = Y)
  
}


.DSA.benchmark.glm <- function(dta, times) {
  fx <- function () {
    for (i in 1:times)
      res <- .DSA.do.lr(dta$x, dta$y)
  }

  gx <- function () {
    for (i in 1:times)
      res <- glm.fit(dta$x, dta$y, family=binomial())
  }
   
  ftran <- system.time(fx())
  gc()
  r <- system.time(gx())

  paste("R:", r[3], "ftran:", ftran[3])
  
}

.DSA.benchmark.lm <- function(dta, times) {
  fx <- function () {
    for (i in 1:times)
      res <- .DSA.do.ols(dta$x, dta$y)
  }

  gx <- function () {
    for (i in 1:times)
      res <- lm.fit(dta$x, dta$y)
  }

  
  ftran <- system.time(fx())
  gc()
  r <- system.time(gx())

  paste("R:", r[3], "ftran:", ftran[3])
}



.DSA.compare.with.R <- function(n = 1000, tol = 1e-10, outcome = 'binary', weighted = FALSE, 
                                their.method = glm.fit, my.method = .DSA.do.lr) {
  xy <- .DSA.generate.X(n, outcome)
  weights <- runif(n)

  if (outcome == 'binary') {
    family = binomial()
  }
  else {
    family = gaussian()
  }
  
  theirs <- c(quote(their.method(xy$x, xy$y, family = family)$dev),
              quote(their.method(xy$x[,1:3], xy$y, family = family)$dev))

  
  ours <- c(quote(my.method(xy$x, xy$y)$dev),
            quote(my.method(xy$x, xy$y, p = 3)$dev))


  if (weighted) {
    theirs <- c(theirs,              
                quote(their.method(xy$x[,1:3], xy$y, weights = weights < .67, family = family)$dev),
                quote(their.method(xy$x, xy$y, family = family, weights = weights < .67)$dev),
                quote(their.method(xy$x, xy$y, family = family, weights = weights)$dev))

    ours <- c(ours,
              quote(my.method(xy$x, xy$y, p = 3, weights = weights < .67)$dev),
              quote(my.method(xy$x, xy$y, weights = weights < .67)$dev),
              quote(my.method(xy$x, xy$y, weights = weights)$dev))
  }


  exprs <- cbind(theirs, ours)

  apply(exprs, MARGIN = 1, function(row) {
    x = sum((eval(row[[1]]) - eval(row[[2]])) > tol)

    if (x != 0) {
      st <- paste("", row[[2]], collapse = "", sep = "")
      print(paste("encountered mismatching values!: ", st))
      print(eval(row[[1]]))
      print(eval(row[[2]]))
    }
    x
  })
  
}
