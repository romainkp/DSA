candidateReduction <- function(formula, data, id = 1:nrow(data), family = gaussian, weights=NULL,...)
{
  extra.args <- list(...)
  if(is.null(extra.args$internal.call) && !is.null(tryCatch(terms(formula),error=function(...){})))
    warning("Except for the outcome, all other variables specified in formula are ignored. To avoid this message in the future, formula should be of the type 'Y ~ .'.\n\n")
  model.mats <- get.model.matrices(formula, data)
  Xlearn <- model.mats$X
  Ylearn <- model.mats$Y
  Y.orig <- model.mats$Y.orig[,1]
  nvarX <- ncol(Xlearn)
  nlearn <- nrow(Xlearn)

  ## now check that we have all the parameters correct and we are ready to
  ## sanely enter the C routines.
  if (!is.numeric(id) || (length(id) != nrow(Xlearn)))
    stop("\nWARNING: id must be an integer vector with length == nrow(data)")
  if (!is.null(weights))
    {
      if (!is.vector(weights) || !is.numeric(weights))
        stop("\nWARNING: weights must be a vector of real numbers.")
      if(length(weights)!=length(Y.orig)[[1]])
        stop("\nWARNING: the length of weights must match the number of observations.\n")
    }

  ##
  ## this was directly lifted from 'glm'
  if (is.character(family)) 
    family <- get(family, mode = "function", envir = parent.frame())
  if (is.function(family)) 
    family <- family()
  if (is.null(family$family)) {
    print(family)
    stop("\nWARNING: 'family' not recognized.\n")
  }
  if ( (family$family != "gaussian") && (family$family != "binomial") && (family$family != "multinomial"))
    stop("\nWARNING: 'family' not recognized.\n")

  Ylearn2 <- Ylearn
  weights2 <- weights
  if(ncol(Ylearn)==2)
    {
      Ntrials <- apply(Ylearn,1,sum)
      if(is.null(weights))weights2 <- Ntrials
      else weights2 <- weights*Ntrials
      Ylearn2 <- Ylearn[,1,drop=FALSE]/Ntrials
    }
  if(!is.null(weights2))weights2[is.na(weights2)] <- 0

  unordered.pval <- matrix(NA,ncol=2,nrow=nvarX)
  dimnames(unordered.pval) <- list(dimnames(Xlearn)[[2]],list("zval/chival","pval"))
  for(i in 1:nvarX)
      unordered.pval[i,] <- get.pval.univariate.gee.independence(Ylearn2,Y.orig,Xlearn[,i],id,family,weights2)
  ordered.pval <- unordered.pval[order(unordered.pval[,"pval"]),,drop=FALSE]
  return(ordered.pval)
}

crossValidate <- function(formula, data, ...) {
  mc <- match.call(expand.dots = TRUE)
  m <- match(c("formula", "data", "id", "family", "weights",
               "usersplits", "userseed", "vfold", "nsplits", "silent"), names(mc), 0)
  mc <- mc[c(1, m)]
  mc[[1]] <- as.name("DSA")

  if (missing(data)) {
    data <- model.frame(formula, environment(formula))
  }
  
  ## get the correct size using the design matrix. 
  mm <- model.matrix(terms(formula, data = data), data = data)[, -1, drop=FALSE]
  size <- length(colnames(mm))
  mc$maxsize <- max(size,1)

  mc$maxsumofpow <- mc$maxorderint <- 1 #not used so okay to set them to 1 now
  mc$cross.validate.model.selected <- FALSE
  mc$expand.forced.term.factors <- TRUE
  mc$data <- data
  if(is.null(mc$silent))mc$silent <- TRUE
  
  res <- eval(mc, parent.frame())
  if(!mc$silent)
    {
      cat("\nThe CV risk for:")
      if(size!=0)print(res$model.selected)
      else print(res$models.allsizes[[1]])
      cat("\nis :",res$average.CVrisks[1, size + 1],".\n")
    }     
  return(res$average.CVrisks[1, size + 1])
}


DSA <- function(formula, data, id = 1:nrow(data), family = gaussian, weights = NULL,candidate.rank = NULL,
                rank.cutoffs =  NULL, maxsize, maxorderint, maxsumofpow, Dmove=TRUE, Smove=TRUE,
                usersplits = NULL,userseed = NULL, vfold = 5, nsplits = 1, model.matrices = FALSE,
                silent = TRUE, ...)
{
  start.time <- as.POSIXlt(Sys.time())
  extra.args <- list(...)
  call <- match.call()
  tryCatch(terms(formula),error=function(...){stop("\n\nWARNING: incorrect formula. The base model must contain the intercept.\n\n")})
  orig.maxsize <- maxsize
  temp <- rep(1,1000)
  storage.mode(temp) <- "integer"
  first.primes <- .C("DSA_PACK_getPrimes",temp)[[1]]
  
  model.mats <- get.model.matrices(formula, data)
  all.vars <- colnames(model.mats$X)

  ## this was directly lifted from 'glm'
  if (is.character(family)) 
    family <- get(family, mode = "function", envir = parent.frame())
  if (is.function(family)) 
    family <- family()
  if (is.null(family$family)) {
    print(family)
    stop("'family' not recognized")
  }

  rank.cutoffs.orig <- rank.cutoffs
  #NO dimension reduction
  if(is.null(rank.cutoffs)) {
    rank.cutoffs <- 1
    candidate.rank <- seq(0,1,by=1/length(all.vars))[-1]
    names(candidate.rank) <- all.vars
  }
  else {   #Dimension reduction
    if(is.null(candidate.rank)) {
      candidate.rank.all <- candidateReduction(formula = formula, data = data, id = id, family = family, weights = weights[nrow(weights),],internal.call=TRUE)
      candidate.rank <- candidate.rank.all[,"pval"]
      if(length(candidate.rank)==1)names(candidate.rank) <- rownames(candidate.rank.all)
    }
    else {
      if(!is.vector(candidate.rank) || !is.numeric(candidate.rank))
        stop("\nWARNING: candidate.rank must be a vector of real values.\n")
      if(any(is.na(candidate.rank)))
        stop("\nWARNING: missing values not allowed in candidate.rank.\n")          
      if(is.null(names(candidate.rank)))
        stop("\nWARNING: the elements of candidate.rank must be named. Their names must correspond to the names of the candidate variables supplied in data.\n")

      #so that ranks are also assigned to the expanded dummy variables (factors)
      if(sum(unlist(lapply(data,is.factor)))!=0) { 
        complete.formula <- as.formula(paste(deparse(as.list(formula)[[2]])," ~ .", sep = ""))
        complete.mf <- model.frame(complete.formula, data = data,na.action = na.pass)
        mt <- attr(complete.mf, "terms")
        for(i in names(attr(mt,"dataClasses")[attr(mt,"dataClasses")=="factor"])) {
          newnames <- colnames(model.matrix(mt[(1:length(attr(mt,"term.labels")))[attr(mt,"term.labels")==i]],
                                            complete.mf, contrasts.arg =
                                            get.contrast.list(as.formula(paste(deparse(as.list(formula)[[2]])," ~ ", i,sep = "")), data))[,-1,drop=FALSE])
          if(is.na(candidate.rank[i]))
            stop("\nWARNING: Missing or unrecognized rank for a factor variable in data.\n")
          newranks <- rep(candidate.rank[i],length(newnames))
          names(newranks) <- newnames
          candidate.rank <- c(candidate.rank[names(candidate.rank)!=i],newranks)
        }
      }
      if(any(!sapply(names(candidate.rank),function(x){sum(all.vars==x)})))
        stop("\nWARNING: Unrecognized variable name in candidate.rank.\n")
      if(length(all.vars)!=length(candidate.rank))
        stop("\nWARNING:  missing rank for at least one candidate variable in data: candidate.rank is nonconformable.\n")
    }
  }
  
  if(!is.vector(rank.cutoffs) || !is.numeric(rank.cutoffs))
    stop("\nWARNING: rank.cutoffs must be a vector of real values.\n")
  if(any(is.na(rank.cutoffs)))
    stop("\nWARNING: missing values not allowed in rank.cutoffs.\n")

  rank.cutoffs <- sort(rank.cutoffs)
  ncutoffs <- length(rank.cutoffs)
  selected.vars <- vector("list",ncutoffs)
  validcutoffs <- NULL #contain the index of valid cut offs stored in rank.cutoffs
  nvalidcutoffs <- 0 #number of valid cut offs
  for(cutoff.loop in 1:ncutoffs) {
    selected.vars[[cutoff.loop]] <- names(candidate.rank[candidate.rank <= rank.cutoffs[cutoff.loop]])
    if(length(selected.vars[[cutoff.loop]])==0)
      warning(paste("\nnone of the candidate variables meets the dimension reduction level corresponding with a rank cut-off of ",
                    rank.cutoffs[cutoff.loop]," - rank cut-off ignored.\n",sep=""))
    else {
      nvalidcutoffs <- nvalidcutoffs+1
      validcutoffs <- c(validcutoffs,cutoff.loop)
    }
  }
  if(nvalidcutoffs==0)
    stop("\nNone of the candidate variables meets the dimension reduction level(s) requested.\n")
  
  if (!missing(silent)) {
    ## get the old one so we can reset it when we leave.
    old.message.level <- getDSAMessageLevel()
    if (silent)
      setDSAMessageLevel(-1)
    else
      setDSAMessageLevel(0)
  }

  ## in this case, we have a model where the factors and the data frame line up;
  ## ie. a call from the DSA to itself in order to cross-validate the final model.
  if (is.null(expand.forced.term.factors <- extra.args$expand.forced.term.factors))
    expand.forced.term.factors <- TRUE
  
  ## in this case, either the DSA is calling itself or someone called the crossValidate
  ## method.
  if (is.null(cross.validate.model.selected <- extra.args$cross.validate.model.selected))
    cross.validate.model.selected <- TRUE

  ## in this case, the DSA is called from DSAgenerate
  if (is.null(CV.DSA <- extra.args$CV.DSA))
    CV.DSA <- TRUE
  
  if (!is.logical(expand.forced.term.factors) ||
      !is.logical(cross.validate.model.selected))
    stop("Error calling DSA, logical parameters must be logical.")
  
  Ylearn <- model.mats$Y
  Xlearn <- model.mats$X
 
  ## save the names from the original.
  Ylearn.orig <- model.mats$Y.orig
  Ylearn.orig.names <- model.mats$Y.names

  ## variables always in the model
  forcedterms <- get.forced.terms(formula, data, Xlearn, expand.forced.term.factors)
  if (!is.null(forcedterms)) {
    nforced <- nrow(forcedterms)
    expanded.formula <- get.final.model(forcedterms,colnames(Xlearn),Ylearn.orig.names) #this is the new formula with factors expanded
  }
  else {
    nforced <- 0
    expanded.formula <- formula
  }
  always.in.Xlearn <- NULL
  if(nforced!=0)
    for(vari in colnames(forcedterms))
      if(sum(forcedterms[,vari])!=0)
        always.in.Xlearn <- c(always.in.Xlearn,vari)

  if(!is.null(always.in.Xlearn))
    always.in.Xlearn <- Xlearn[,always.in.Xlearn,drop=FALSE]

  if (model.matrices == TRUE) {
    return(model.mats)
  }  
  nlearn <- nrow(Xlearn)
    
  if(!is.null(usersplits))
    if(nsplits != 1)nsplits <- 1
  
  ## now check that we have all the parameters correct and we are ready to
  ## sanely enter the C routines.
  if (!is.numeric(id) || (length(id) != nrow(Xlearn)))
    stop("\nWARNING: id must be an integer vector with length == nrow(data)")
  if (!(is.numeric(vfold) && (length(vfold) == 1) && (vfold > 1)))
    stop("\nWARNING: vfold must be a positive integer > 1")
  if (is.null(usersplits) && !is.null(userseed)) {
    if (!(is.numeric(userseed) && length(userseed) == 1))
      stop("\nWARNING: userseed must be a scalar integer")
  }
  if (!(is.numeric(maxorderint) && maxorderint > 0) ||
      !(is.numeric(maxsumofpow) && maxsumofpow > 0))
    stop("\nWARNING: maxorderint and maxsumofpow must be integers strictly greater than 0")
  if(!(is.numeric(maxsize) && maxsize >= 0))
    stop("\nWARNING: maxsize must be an integer greater than 0")
  if(maxsize< nforced)
    stop(paste("\nWARNING: maxsize (=",maxsize,") must be larger than or equal to the number of forced terms (=",nforced,").",sep=""))
  if(maxorderint>maxsumofpow)stop("\nWARNING: maxorderint must be lower or equal to maxsumofpow.\n")
  
  if (!is.null(usersplits) && !(is.matrix(usersplits) && is.numeric(usersplits) &&
                                all(!is.na(usersplits)) && all(usersplits == 1 || usersplits == 0)))
    stop("\nWARNING: usersplits must be a matrix of 0s or 1s.\n")

  if(!is.null(usersplits)) {
    vfold <- dim(usersplits)[[1]] #allows for alternate splitting scheme (not v fold)
    if(dim(usersplits)[[2]] != nrow(Xlearn))
      stop("\nWARNING: the number of columns of usersplits does not correspond to the number of observations.\n")
    if(!all(apply(usersplits,1,sum) > 0 & apply(usersplits,1,sum) < nrow(Xlearn)))
      stop("\nWARNING: the number of observations in each validation sample should be strictly larger than 0 and strictly lower than the total number of observations.")
    }
  if (!is.null(weights))
    {
      if (!is.matrix(weights) || !is.numeric(weights))
        stop("\nWARNING: weights must be a matrix of real numbers.")
      if(dim(weights)[[2]] != nrow(Xlearn))
        stop("\nWARNING: the number of columns in weights must match the number of observations.\n")
      if(is.null(usersplits))
        {
          if(dim(weights)[[1]]!=(vfold+1))
            stop("\nWARNING: the number of rows in weights must be equal to the number of data splits plus one.\n")
        }
      else
        {
          if(dim(weights)[[1]]!=(nrow(usersplits)+1))
            stop("\nWARNING: the number of rows in weights must be equal to the number of data splits plus one.\n")
        }      
    }

  ## consistency between the arguments. (this check only makes sense when the default vfold splitting scheme is used)
  if(is.null(usersplits) && (length(unique(id))/vfold) <= 1)
    stop("\nWARNING: the number of independent experimental units is too low for vfold.")


  ## set the multinomial dimensions (M) for Ylearn.
  mLevels <- if (!is.null(model.mats$class.labels)) {
    length(model.mats$class.labels) - 1
  }
  else {
    1
  }
  
  CVrisks <- matrix(0, nrow = maxorderint, ncol = maxsize+1)
  coefficients <- rep(0, (maxsize + 1) * mLevels)

  if (is.null(weights)) {
    WTtrainvallearn <- matrix(1, ncol = nlearn, nrow = ifelse(is.null(usersplits),vfold+1,nrow(usersplits)+1))
    useWeights <- 0
  }
  else {
    WTtrainvallearn <- weights
    useWeights <- 1
  }
  
  if (ncol(Ylearn) == 2) {
    binWTtrainvallearn <- Ylearn[,1] + Ylearn[,2]
    
    # now Ylearn becomes the success ratio.
    tempYlearn <- as.matrix(ifelse(binWTtrainvallearn == 0, 0, Ylearn[,1]/binWTtrainvallearn))
    dimnames(tempYlearn)[[2]] <- list(dimnames(Ylearn)[[2]][1])
    Ylearn <- tempYlearn
    
    # and the weights become the number of trials:
    binWTtrainvallearn[is.na(binWTtrainvallearn)] <- 0 #first, so that WTtrainvallearn=NA if user weights=NA only

    WTtrainvallearn <- t(apply(WTtrainvallearn, 1, function(row) {
      row * binWTtrainvallearn
    }))

    useWeights <- 1
  }
  else {
    binWTtrainvallearn <- rep(1, nlearn)
  }

  
  binind <- switch(family$family,
                   "multinomial" = mLevels,
                   "gaussian" = 0,
                   "binomial" = 1)

  if(sum(is.na(WTtrainvallearn))>0)stop("\nWARNING: weights cannot be missing.\n")
  
  storage.mode(Ylearn) <- "double"
  storage.mode(coefficients) <- "double"
  storage.mode(CVrisks) <- "double"
  storage.mode(WTtrainvallearn) <- "double"
  storage.mode(binWTtrainvallearn) <- "double"

  vfoldscheme <- FALSE  
  if(is.null(usersplits)) {
    # compute splits. 
    if (!missing(userseed) & is.numeric(userseed))
      set.seed(userseed)
    vfoldscheme <- TRUE
  }

  averageCVs <- array(0,c(maxorderint,maxsize+1,nvalidcutoffs))
  allsplits <- vector("list",nsplits)
  for(isplit in 1:nsplits)
    {
      if(!silent && vfoldscheme && cross.validate.model.selected)cat("\nV-fold split",isplit,"out of",nsplits,"\n")
      if(vfoldscheme) {  
        uniqueIDlearn <- t(as.matrix(unique(id)))
        nIDlearn <- ncol(uniqueIDlearn)
        
        #each unique ID is assigned to a split:
        splitvector <- c(rep(1:vfold, each = floor(nIDlearn/vfold)),
                         rep(vfold, nIDlearn - floor(nIDlearn/vfold)*vfold))
        splitvector <- sample(splitvector)

        #each observation is assigned to the split number to which its
        #corresponding ID was assigned
        splitvectorlearn <- apply(t(id), 2, function(x,y,z) {
          return(y[z==x])
        }, y=splitvector,z=uniqueIDlearn)
        
        usersplits <- apply(t(splitvectorlearn),2,function(x,y) {
          splitassignment <- rep(0,y)
          splitassignment[x] <- 1
          return(splitassignment)
        }, y=vfold) 
      }
      maxnval <- max(apply(usersplits, 1, sum))
      maxntrain <- nlearn - min(apply(usersplits,1,sum))
      storage.mode(usersplits) <- "integer"
      storage.mode(maxnval) <- "integer"    
      storage.mode(maxntrain) <- "integer"    
      storage.mode(nforced) <- "integer"
      
      for(validcutoffs.loop in 1:nvalidcutoffs)
        {
          cutoff.loop <- validcutoffs[validcutoffs.loop]
          Xlearn <- model.mats$X[,selected.vars[[cutoff.loop]],drop=FALSE]
          
          if(!is.null(always.in.Xlearn))
            {
              for(vari in colnames(always.in.Xlearn))if(sum(colnames(Xlearn)==vari)==0)
                {
                  temp.exp <- parse(text=paste('Xlearn <- cbind(Xlearn,"',vari,'"=always.in.Xlearn[,vari])',sep=""))
                  eval(temp.exp)
                }
            }

          nvarX <- ncol(Xlearn)
          ## construct the forced term matrix - null if no forced terms.
          forcedterms <- get.forced.terms(formula, data, Xlearn, expand.forced.term.factors)
          bestmodelsspace <- rep(NA,maxsize*maxsize*nvarX)

          if((nvarX<=1000) && (first.primes[nvarX]^maxsumofpow)<=(2^32-1) )usetree <- 1
          else
            {
              usetree <- 0
              cat("\nWARNING: Due to the large number of candidate variables considered, the DSA algorithm invoked cannot avoid refitting models already evaluated. You may want to reduce the number of candidate variables or further constrain the complexity of the models considered.\n")
            }
          gCONST <- c(maxsize, maxorderint, maxsumofpow, vfold, nlearn,
                      nvarX, useWeights, NA, NA, binind, as.numeric(CV.DSA),
                      usetree,nsplits,as.numeric(Dmove),as.numeric(Smove))
          
          storage.mode(gCONST) <- "integer"
          storage.mode(Xlearn) <- "double"        
          storage.mode(bestmodelsspace) <- "integer"
          storage.mode(forcedterms) <- "integer"

          on.exit(.C("DSA_PACK_DeleteTree"))
          res <- .C("DSA_PACK_Rentry",
                    gCONST,
                    Xlearn,
                    Ylearn,
                    CVrisks,
                    coefficients,
                    usersplits,
                    WTtrainvallearn,
                    nforced,
                    forcedterms,
                    maxntrain,
                    maxnval,
                    bestmodelsspace,
                    binWTtrainvallearn,
                    NAOK = TRUE)
          averageCVs[,,validcutoffs.loop] <- averageCVs[,,validcutoffs.loop]+res[[4]]
          gCONST <- res[[1]]
          if(!silent && vfoldscheme && cross.validate.model.selected)
            {
              cat("\nAverage CV risks for split",isplit,"and cutoff",rank.cutoffs[cutoff.loop],":\n")
              print(res[[4]])
              cat("\nSmallest CV risks for split",isplit,"and cutoff",rank.cutoffs[cutoff.loop],":",min(res[[4]]),"\n")
            }
        }
      allsplits[[isplit]] <- usersplits
    }



  if(nsplits==1 && nvalidcutoffs==1)
    {
      models.allsizes <- vector("list",maxsize+1)
      if(nforced==0)
        models.allsizes[[1]] <- get.final.model(matrix(0,nrow=maxsize,ncol=nvarX),
                                                dimnames(Xlearn)[[2]], Ylearn.orig.names)
      for(i in max(1,nforced):maxsize) {
        model.sizei <- matrix(res[[12]][(1+(i-1)*(maxsize*nvarX)):(i*(maxsize*nvarX))],ncol=nvarX)
        if(sum(is.na(model.sizei))==0)
          models.allsizes[[i+1]] <- get.final.model(model.sizei, dimnames(Xlearn)[[2]], Ylearn.orig.names)
      }
      modelf <- models.allsizes[[gCONST[8]+1]]
      coefficients <- res[[5]]
      bestcutoff <- 1
    }
  else
    {
      averageCVs <- averageCVs/nsplits
      minpos <- which(averageCVs==min(averageCVs),arr.ind=TRUE)
      if(nrow(minpos)>1  && cross.validate.model.selected)warning("\nThere were ties in the array of crossvalidated risks.")
      minpos <- minpos[order(minpos[,"dim3"],minpos[,"dim1"]),,drop=FALSE]
      bestorderint <- minpos[1,"dim1"]
      bestsize <- minpos[1,"dim2"]
      bestcutoff <- minpos[1,"dim3"]

      gCONST[8] <- bestsize-1
      gCONST[9] <- bestorderint
      
      if(!silent && cross.validate.model.selected)
        cat("\nModel selected of size ",bestsize," and order of interaction ",
            bestorderint," based on the rank cut-off ",rank.cutoffs[validcutoffs[bestcutoff]],".\n",sep="")

      bestXlearn <- model.mats$X[,selected.vars[[validcutoffs[bestcutoff]]],drop=FALSE]
      if(!is.null(always.in.Xlearn))
        {
          for(vari in colnames(always.in.Xlearn))if(sum(colnames(bestXlearn)==vari)==0)
            {
              temp.exp <- parse(text=paste('bestXlearn <- cbind(bestXlearn,"',vari,'"=always.in.Xlearn[,vari])',sep=""))
              eval(temp.exp)
            }
        }
      bestdata <- as.data.frame(cbind(Ylearn.orig,bestXlearn))
      finalDSA <- DSAgenerate(formula=expanded.formula, data = bestdata, family = family, weights = weights,
                              maxsize = maxsize, orderint = bestorderint, maxsumofpow = maxsumofpow,
                              silent = silent, expand.forced.term.factors = FALSE)
      models.allsizes <- finalDSA$models.allsizes
      modelf <- models.allsizes[[bestsize]]
      mtofit <- list(family=finalDSA$family,formula = formula,data = bestdata, weights = weights, model.selected = modelf)
      coefficients <- as.vector(t(coefficients(get.dsa.fit(mtofit))))
    }
  

  ## now we will cross validate the response. calling the DSA again w/appropriate
  ## flags set. This is similar to 'crossValidate' however, the manner in which the
  ## size is obtained is a bit simpler because the DSA only returns numeric columns.
  cv.res <- NA #need to keep this initialization in case DSA is called with cross.validate.model.selected=FALSE
  if (cross.validate.model.selected) {
    if(!silent)cat("\nComputing the CV risk for the final model\n")
    size <- length(labels(terms(modelf, data = data)))
    if (size > 0) {
      cv.call <- call
      cv.call$formula <- modelf
      cv.call$maxsize <- size
      cv.call$cross.validate.model.selected <- FALSE
      cv.call$expand.forced.term.factors <- FALSE
##       cv.call$rank.cutoffs <- rank.cutoffs[validcutoffs[bestcutoff]]
      cv.call$rank.cutoffs <- NULL

      cv.res <- try(eval(cv.call, parent.frame()),silent=TRUE)
      if(class(cv.res)=="try-error")
        {
          cat(cv.res)
          cat("\nThe CV risk for:")
          print(modelf)
          cat("\ncould not be computed for the data split considered.\n")
          cv.res <- NA
        }
      else
        {
          if(!silent)
            {
              cat("\nThe CV risk for:")
              print(cv.res$model.selected)
              cat("\nis :",cv.res$average.CVrisks[1, size + 1],".\n",sep="")
            }
          cv.res <- cv.res$average.CVrisks[1, size + 1]
        }
    }
    else { #if the interecpt model is selected we can read off the CV risk from averageCVs
      cv.res <- averageCVs[1,1,1]
      if(!silent)
        {
          cat("\nThe CV risk for:")
          print(models.allsizes[[1]])
          cat("\nis :",cv.res,".\n")
        }
    }
  }

  if(nsplits==1)usersplits <- allsplits[[1]]
  else usersplits <- allsplits
  
  candidate.vars <- selected.vars[[cutoff.loop]]
  if(!is.null(always.in.Xlearn)) #append forced variables
    {
      for(vari in colnames(always.in.Xlearn))if(sum(candidate.vars==vari)==0)
        {
          temp.exp <- parse(text=paste('candidate.vars <- c(candidate.vars,"',vari,'")',sep=""))
          eval(temp.exp)
        }
    }
  if(dim(averageCVs)[3]==1)
    {
      averageCVs <- matrix(as.vector(averageCVs[,,,drop=TRUE]),nrow = maxorderint, ncol = maxsize+1)
      colnames(averageCVs) <- paste("size",1:(maxsize+1))
      rownames(averageCVs) <- paste("interaction",1:maxorderint)
      min.risk <- c(gCONST[9], gCONST[8]+1)
      names(min.risk) <- c("interaction","size")
      selected.vars <- candidate.vars
    }
  else
    {
      colnames(averageCVs) <- paste("size",1:(maxsize+1))
      rownames(averageCVs) <- paste("interaction",1:maxorderint)
      dimnames(averageCVs)[[3]] <- paste("rank cut-off",rank.cutoffs[validcutoffs])
      min.risk <- c(gCONST[9], gCONST[8]+1, bestcutoff)
      names(min.risk) <- c("interaction","size","cut-off")
      selected.vars <- selected.vars[[validcutoffs[bestcutoff]]]
      if(!is.null(always.in.Xlearn)) #append forced variables
        {
          for(vari in colnames(always.in.Xlearn))if(sum(selected.vars==vari)==0)
            {
              temp.exp <- parse(text=paste('selected.vars <- c(selected.vars,"',vari,'")',sep=""))
              eval(temp.exp)
            }
        }
    }
  if(is.null(rank.cutoffs.orig) || (nvalidcutoffs==1))selected.vars <- NULL #no selection if no dimension reduction
  coef.matrix <- matrix(coefficients[1:((gCONST[8] + 1)*mLevels)], nrow = mLevels, byrow = TRUE)
  rownames(coef.matrix) <- model.mats$class.labels[-1] ## chop off the baseline category.
  if(!is.null(modelf))colnames(coef.matrix) <- c("(Intercept)",attributes(terms(modelf))$term.labels)

  final.time <- as.POSIXlt(Sys.time())
  total.time <- difftime(final.time, start.time)
  
  output <- list("average.CVrisks" = averageCVs, "min.risk" = min.risk, "final.cutoff"=rank.cutoffs[validcutoffs[bestcutoff]],
                 "model.selected" = modelf, "family" = family$family,
                 "coefficients" = coef.matrix,
                 "average.CVrisk.model.selected" = cv.res, "models.allsizes"=models.allsizes,
                 "candidate.vars"=candidate.vars, "selected.vars"=selected.vars,
                 "candidate.rank"=candidate.rank, #number of variable passed in to the DSA
                 "rank.cutoffs"=rank.cutoffs.orig, "computing.time" = total.time, "call" = call,
                 "formula" = formula,"data" = data, "id" = id, "userseed" = userseed,
                 "splits" = usersplits, "weights" = weights, "moves" = c("Dmove"=Dmove,"Smove"=Smove) )
  
  class(output) <- "DSA"

  ## reset some things on the way out.
  if (!missing(silent)) 
    setDSAMessageLevel(old.message.level)
  return(output)
}

DSAgenerate <- function(formula, data, family = gaussian, weights = NULL, maxsize, orderint,
                        maxsumofpow, silent = TRUE, ...)
{
  extra.args <- list(...)    
  call <- match.call()
  if(!is.null(weights))
    {
      if(!is.matrix(weights))weights <- matrix(weights,nrow=1)
      DSAweights <- matrix(weights[nrow(weights),],nrow=6,ncol=ncol(weights),byrow=T)
    }
  else DSAweights <- NULL
  if (is.null(expand.forced.term.factors <- extra.args$expand.forced.term.factors))
    expand.forced.term.factors <- TRUE  

  res <- DSA(formula=formula, data=data, family=family, weights=DSAweights,
             maxsize=maxsize, maxorderint=orderint, maxsumofpow=maxsumofpow, silent=silent,
             CV.DSA = FALSE, cross.validate.model.selected = FALSE,
             expand.forced.term.factors = expand.forced.term.factors)
  
  output <- list("models.allsizes"=res$models.allsizes, "family" = res$family,
                 "computing.time" = res$computing.time, "call" = call)
  return(output)
}

setDSAMessageLevel <- function(newlevel) {
  storage.mode(newlevel) <- "integer"
  res <- .C("DSA_PACK_setMessagelevel", newlevel)[[1]]
}

getDSAMessageLevel <- function() {
  .C("DSA_PACK_getMessageLevel", integer(1))[[1]]
}

getDSASubversionInfo <- function() {
  svn.file <- paste(system.file(package = "DSA"), "/", "svn.info", sep = "")
  if (!file.exists(svn.file))
    stop("Package built without subversion information.")
  else
    xx <- readLines(svn.file)

  print.package.message(xx)
  
}

print.package.message <- function(svn.info = NULL) {
  if (!is.null(svn.info))
    svn.info <- strsplit(svn.info, ": ")

  cat("\n DSA version: ")
  cat(packageDescription('DSA')$Version)
  cat("\n Subversion Repository Version: ")
  if (!is.null(svn.info)) {
    cat(strsplit(svn.info[[1]][9], " ")[[1]][1])
    cat("\n Last changed date: ")
    cat(svn.info[[1]][10])
  }
  else {
    cat("Development Version")
  }
  cat("\n\n")
}

hasDSASubversionInfo <- function() {
  svn.file <- paste(system.file(package = "DSA"), "/", "svn.info", sep = "")
  if (!file.exists(svn.file))
    FALSE
  else
    TRUE
}

####################################################################################
### ################################ S3 Methods ################################ ###
####################################################################################

print.DSA <- function(x, ...) {
  cat("\nModel Selected:\n\t")
  print(x$model.selected)
##   cat("\nModel Coefficients:\n")
##   print(x$coefficients)
  cat("\nTime:\n\t", x$computing.time,attributes(x$computing.time)[[2]])
  cat("\n\n")
}

predict.DSA <- function(object, newdata = NULL, ...) {
  fit <- get.dsa.fit(object)

  if (missing(newdata)) {
    predict(fit, ...)
  }
  else {
    predict(fit, newdata, ...)
  }
}

plot.DSA <- function(x, plot.compare = FALSE, ...) {
  compare.plot <- function(dsa.res, ...) {
    nplots <- prod(dim(dsa.res$average.CVrisks)[-c(1,2)])
    par(mfrow=n2mfrow(nplots))
    ylimits <- range(dsa.res$average.CVrisks[which(is.finite(dsa.res$average.CVrisks))])
    if(nplots==1)
      {
        if(nrow(dsa.res$average.CVrisks)>1)
          {
            plot(dsa.res$average.CVrisks[1,], xlab = "Model Size",
                 ylab = "Average CV Risks",
                 pch = 1,ylim=ylimits)
            
            for (i in 2:nrow(dsa.res$average.CVrisks)) {
              points(dsa.res$average.CVrisks[i,], pch = i)
            }
            avg.risks <- dsa.res$average.CVrisks[which(is.finite(dsa.res$average.CVrisks))]
            avg.risks <- range(avg.risks)
            avg.risks <- (3*avg.risks[2]+avg.risks[1])/4
            ##         avg.risks <- avg.risks[2] - (avg.risks[2] - avg.risks[1])*.80 
            legend(ncol(dsa.res$average.CVrisks)*.6, avg.risks,
                   pch = 1:nrow(dsa.res$average.CVrisks),
                   legend = paste("interaction order", 1:nrow(dsa.res$average.CVrisks)))
          }
        else
          {
            plot.DSA(dsa.res, plot.compare= FALSE,...)
          }        
      }
    else
      {
        for(j in 1:nplots)
          {
            plot(dsa.res$average.CVrisks[1,,j], xlab = "Model Size",
                 ylab = "Average CV Risks", main = paste(dimnames(x$average.CVrisks)[[3]][j]),
                 pch = 1,ylim=ylimits)
            
            for (i in 2:nrow(dsa.res$average.CVrisks)) {
              points(dsa.res$average.CVrisks[i,,j], pch = i)
            }
            avg.risks <- dsa.res$average.CVrisks[which(is.finite(dsa.res$average.CVrisks))]
            avg.risks <- range(avg.risks)
            avg.risks <- (3*avg.risks[2]+avg.risks[1])/4
            ##             avg.risks <- avg.risks[2] - (avg.risks[2] - avg.risks[1])*.80 
            legend(ncol(dsa.res$average.CVrisks)*.6, avg.risks,
                   pch = 1:nrow(dsa.res$average.CVrisks),
                   legend = paste("interaction order", 1:nrow(dsa.res$average.CVrisks)))
          }
      }
  }
  
  if (plot.compare)
    compare.plot(x, ...)
  else {
    orig.mfrow <- par()$mfrow
    nplots <- prod(dim(x$average.CVrisks)[-2])
    par(mfrow=n2mfrow(nplots))
    if(length(dim(x$average.CVrisks))==2)
      {
        for (i in 1:nrow(x$average.CVrisks)) {
          plot(x$average.CVrisks[i,], xlab = "Model Size",
               ylab = "Average CV Risks", main = paste("Interaction order",i))
        }
      }
    else
      {
        for(j in 1:dim(x$average.CVrisks)[[3]])for (i in 1:nrow(x$average.CVrisks)) {
          plot(x$average.CVrisks[i,,j], xlab = "Model Size",
               ylab = "Average CV Risks", main = paste("Interaction order",i,"and",dimnames(x$average.CVrisks)[[3]][j]))
        }
      }
    par(mfrow=orig.mfrow)
  }
}

summary.DSA <- function(object, ...) {
  fit <- get.dsa.fit(object)
  coefficients <- coefficients(fit)
  if(is.matrix(coefficients))colnames(coefficients) <- c("(Intercept)",attributes(terms(object$model.selected))$term.labels)
  
  mm <- get.model.matrices(object$formula, data = object$data)
  df <- as.data.frame(cbind(mm$Y.orig, mm$X))
  n.obs <- nrow(df)
  n.obs.nonna <- nrow(model.matrix(object$model.selected, df))
  
  n.id <- length(unique(object$id))
  nid.mf <- model.frame(object$model.selected, data = df, na.action = na.pass)
  nid.mt <- attr(nid.mf, "terms")
  nid.mm <- model.matrix(nid.mt, nid.mf)[, -1]
  nid.mr <- model.response(nid.mf)
  mat <- cbind(nid.mr, nid.mm)
#  mat <- cbind(id = object$id, nid.mr, nid.mm)

  non.na.obs <- apply(mat,1,function(x){return(!any(is.na(x)))})
  non.na.obs.per.id <- unlist(lapply(split(non.na.obs,object$id),sum))
  n.id.nonna <- sum(non.na.obs.per.id!=0)
    
##   n.id.nonna <- 0
##   for (mid in unique(object$id)) {
##     if (!all(is.na(rowSums(mat[mat[,1] == mid, , drop = FALSE]))))
##       n.id.nonna <- n.id.nonna + 1
##   }

  if(length(dim(object$average.CVrisks))==2)
    {
      cvrisk <- object$average.CVrisks[object$min.risk[1], object$min.risk[2]]
      best.rank.cutoff <- NULL
    }
  else
    {
      best.rank.cutoff <- object$final.cutoff
      cvrisk <- object$average.CVrisks[object$min.risk[1], object$min.risk[2],object$min.risk[3]]
    }
  io.and.size <- list("model.size.selected" = object$min.risk[2],
                      "interaction.order.selected" = object$min.risk[1])

  summary.DSA <- list(model.selected = object$model.selected,
                      family = object$family,
                      n.obs = n.obs,
                      n.obs.nonna = n.obs.nonna,
                      best.average.CVrisk = cvrisk,
                      coefficients = coefficients,
                      io.and.size = io.and.size,
                      average.CVrisk.model.selected = object$average.CVrisk.model.selected,
                      n.id = n.id, n.id.nonna = n.id.nonna, nselected.vars = length(object$selected.vars),
                      ncandidate.vars = length(object$candidate.vars), n.vars = length(object$candidate.rank),
                      rank.cutoffs =object$rank.cutoffs,best.rank.cutoff=best.rank.cutoff,formula=object$formula,
                      moves = object$moves)
  
  class(summary.DSA) <- "summary.DSA"
  return(summary.DSA)
}

print.summary.DSA <- function(x, ...) {
  cat("\nModel selected:\n\t")
  print(x$model)
  cat("\n")
  print(x$coefficients)
  cat("\nBase model: ")
  print(x$formula)
  cat("Family: ")
  cat(x$family)
  if(all(x$moves==c(FALSE,TRUE)))cat("\nDeletion moves inhibited")
  else if(all(x$moves==c(TRUE,FALSE)))cat("\nSubstitution moves inhibited")
  else if(all(x$moves==c(FALSE,FALSE)))cat("\nDeletion and substitution moves inhibited")  
  cat("\nTotal number of experimental units: ")
  cat(x$n.id)
  cat("\nNumber of experimental units used to fit the model selected: ")
  cat(x$n.id.nonna)
  cat("\nTotal number of observations: ")
  cat(x$n.obs)
  cat("\nNumber of observations used to fit the model selected: ")
  cat(x$n.obs.nonna)
  
  cat("\n\nDimension reduction: ")  
  if(is.null(x$rank.cutoffs))
    {
      cat("no")
      cat("\nNumber of candidate variables considered:",x$n.vars)
      cat("\nMinimum average cross-validated risk for model size ", x$io.and.size$model.size.selected,
          ", interaction order ", x$io.and.size$interaction.order.selected,": ",sep="")
      cat(x$best.average.CVrisk)      
    }
  else
    {
      cat("yes")
      cat("\nNumber of candidate variables supplied:",x$n.vars)      
      if(is.null(x$best.rank.cutoff))
        {
          cat("\nUser-specified rank cut-off: ",x$rank.cutoffs," (",x$ncandidate.vars," variables considered)",sep="")
          rank.temp <- x$rank.cutoffs
        }
      else
        {
          cat("\nNumber of candidate variables considered:",x$ncandidate.vars)
          cat("\nSelected rank cut-off: ",x$best.rank.cutoff," (",x$nselected.vars," variables)",sep="")
          rank.temp <- x$best.rank.cutoff
        }
      cat("\nMinimum average cross-validated risk for model size ", x$io.and.size$model.size.selected,
          ", interaction order ", x$io.and.size$interaction.order.selected, " and rank cut-off ",rank.temp,": ",sep="")
      cat(x$best.average.CVrisk)      
    }
  
##   if((x$ncandidate.vars==x$n.vars))
##     {
##       cat("none")
##     }
##   cat("\n\tNumber of candidate variables:",x$n.vars)      
##   cat("\n\tNumber of candidate variables truly considered:",x$ncandidate.vars)
##   if(is.null(x$best.rank.cutoff))
##     {
##       cat("\n\tUser-specified rank cut-off: ",x$rank.cutoffs," (",x$nselected.vars," variables)",sep="")
##       rank.temp <- x$rank.cutoffs
##     }
##   else
##     {
##       cat("\n\tSelected rank cut-off: ",x$best.rank.cutoff," (",x$nselected.vars," variables)",sep="")
##       rank.temp <- x$best.rank.cutoff
##     }

  cat("\nAverage cross-validated risk for the model selected: ")
  cat(x$average.CVrisk.model.selected)
  cat("\n\n")
}

coefficients.DSA <- function(object, ...) {
  mm.final <- model.matrix(object$model.selected, data = object$data)
  colnames(object$coefficients) <- colnames(mm.final)
  return(object$coefficients)
}

residuals.DSA <- function(object, ...) {
  fit <- get.dsa.fit(object, ...)
  return(residuals(fit, ...))
}

#######################################################################################
### ############################### non-public methods ############################ ###
#######################################################################################
get.pval.univariate.gee.independence <- function(Y,Y.orig,X,id,family,weights)
  {
    if(missing(weights) || is.null(weights))weights <- rep(1,length(X))
    else if(sum(is.na(weights))!=0)stop("WARNINGS: weights cannot be missing.")

    if(family$family == "multinomial")
      {
        Yvec <- matrix(Y,byrow=TRUE,ncol=length(X))
        m <- length(Y)/length(X)

        dat <- cbind(X,t(Yvec))
        obs.nonna <- apply(dat,1,function(x){return(sum(is.na(x))==0)})
        weights[!obs.nonna] <- 0
        X[!obs.nonna] <- 0
        Yvec[,!obs.nonna] <- 0
#       id.nonna <- apply(cbind(unique(id)),1,function(x,y,z){return(sum(y[z==x]))},y=obs.nonna,z=id)
        non.na.obs.per.id <- unlist(lapply(split(obs.nonna,id),sum))
        n <- sum(non.na.obs.per.id!=0)
        
#        n <- sum(id.nonna!=0)

        coef <- get.gee.multinom(Yvec=Yvec,Y.orig=Y.orig,X=X,m=m,id=id,weights=weights,n=n)        
        IC.n <- getIC(Yvec,X,coef,m=m,id=id,weights=weights,n=n)
        chival <- get.test.stat(coef,IC.n,m=m,n=n)
        pval <- pchisq(chival,df=m,lower=FALSE)
        res <- cbind(chival,pval)
      }
    else
      {
        coef <- glm(Y~X,family=family,weights=weights)$coef
#        cat("\nFrom candidate - coef:",coef[2])
        
        dat <- cbind(X,Y)
        deltaxy <- apply(dat,1,function(x){return(sum(is.na(x))==0)})
#       id.nonna <- apply(cbind(unique(id)),1,function(x,y,z){return(sum(y[z==x]))},y=deltaxy,z=id)
        non.na.obs.per.id <- unlist(lapply(split(deltaxy,id),sum))
        nobs <- sum(non.na.obs.per.id!=0)
        
#        nobs <- sum(id.nonna!=0)        
#        uid <- unique(id)[id.nonna!=0]
        if(family$family=="gaussian")
          {
            epsilon.v <- as.vector(Y-cbind(1,X)%*%cbind(coef))
            sqrt.C.v <- deltaxy
          }
        else if(family$family=="binomial")
          {
            temp <- exp(-as.vector(cbind(1,X)%*%cbind(coef)))
            epsilon.v <- as.vector(Y-1/(1+temp))
            sqrt.C.v <- sqrt(temp/(1+temp)^2)
          }          
        X.v <- X
        epsilon.deltaxy <- epsilon.v
        X.deltaxy <- X.v

##         missing.obs <- as.logical(is.na(epsilon.v) + is.na(X.v) + is.na(sqrt.C.deltaxy)) # same as !deltaxy
        epsilon.deltaxy[!deltaxy] <- 0
        X.deltaxy[!deltaxy] <- 0
        weights[!deltaxy] <- 0
        sqrt.C.deltaxy <- sqrt.C.v*sqrt(weights)        
        sqrt.C.deltaxy[!deltaxy] <- 0
        
#        epsilon.deltaxy[is.na(epsilon.v)] <- 0
#        X.deltaxy[is.na(X.v)] <- 0
#        sqrt.C.deltaxy[is.na(sqrt.C.deltaxy)] <- 0
#        EF.v <- rbind(weights*epsilon.deltaxy, weights*X.deltaxy*epsilon.deltaxy)
#        EF <- apply(cbind(uid),1,function(x){return(apply(cbind(EF.v[,id==x]),1,sum))})
        EF.v <- cbind(weights*epsilon.deltaxy, weights*X.deltaxy*epsilon.deltaxy)        
        EF <- t(rowsum(EF.v,id))
          

        temp.mat <- rbind(sqrt.C.deltaxy,X.deltaxy*sqrt.C.deltaxy)
        NormC.inv <- temp.mat%*%t(temp.mat)/nobs
#        NormC <- chol2inv(chol(NormC.inv))
        NormC <- solve(NormC.inv)
        IC <- NormC%*%EF
        var.beta <- IC%*%t(IC)/nobs
        se.beta <- sqrt(diag(var.beta)/nobs)
#        cat("\nFrom candidate - sd:",se.beta[2])
        zval <- coef/se.beta

        pval <- 2 * pnorm(abs(zval), lower.tail = FALSE)
        res <- cbind(zval,pval)[2,]
      }
    
    return(res)
  }

get.gee.multinom <- function(Yvec,Y.orig,X,m,id,weights=NULL,conv.criteria=1e-16,linear.search=FALSE,n,silent=TRUE)
  {
    res <- fitModel(Y.orig ~ X, family=multinomial,baseline=1,weights=weights)
    resk1.test <- cbind(as.vector(t(res$estimates)))
    j <- 0
    srule.n.old <- Inf
    epsilon <- 0
    while(1)
      {
        j <- j+1
        IC.n <- getIC(Yvec,X,resk1.test,m,id,weights,n)
        c.n <- apply(IC.n,1,mean)
        srule.n <- sum(c.n^2)
        if(!silent)
          {
            cat("\tEstimate tried",j,"-",epsilon,":")
            for(i in 1:m)
              {
                cat("\n\t\t",resk1.test[((i-1)*2+1):(i*2)])
              }
            cat(" - stopping rule: ",srule.n,sep="")
          }
        if(srule.n<conv.criteria)break
        else
          {
            if(srule.n>srule.n.old)
              {
                if(!silent)cat(" - rejected\n")                
                if(linear.search)epsilon <- epsilon+0.1
                else epsilon <- 1
                if(abs(epsilon-1)<1e-8 )
                  {
                    cat("\nNo Convergence\n")
                    return(NULL)
                  }
                resk1.test <- epsilon*resk+(1-epsilon)*resk1.try
              }
            else
              {
                if(!silent)cat(" - accepted\n")                
                resk <- resk1.test
                epsilon <- 0
                resk1.test <- resk1.try <- resk + c.n
                srule.n.old <- srule.n
              }
          }
      }
    return(resk1.test)
  }

phi <- function(X,coef.k,m)
{
  denom <- 1
  res <- matrix(NA,nrow=m,ncol=ncol(X))  
  for(i in 1:m)
    {
      res[i,] <- exp(rbind(coef.k[((i-1)*2+1):(i*2),1])%*%X)
      denom <- denom+res[i,]
    }
   res <- sweep(res,2,denom,"/")
  
  return(res)
}

getIC <- function(Yvec,X,coef.k,m,id,weights=NULL,n)
  {
    nobs <- ncol(Yvec)

    #estimating function
    X <- rbind(1,matrix(X,ncol=nobs))
    
    phi.p <- phi(X,coef=coef.k,m)
    eps <- Yvec-phi.p

    A <- apply(X,2,function(x,y){return(rep(x,y))},y=m)
    #weights:
    if(!is.null(weights))
      {
        A.w <- sweep(X,2,weights,"*")
        A <- apply(A.w,2,function(x,y){return(rep(x,y))},y=m)
      }
    B <- apply(phi.p,2,function(x,y){return(rep(x,each=y))},y=2)
    C <- apply(eps,2,function(x,y){return(rep(x,each=y))},y=2)
    M <- A*B*C

    N <- matrix(0,nrow=m*2,ncol=nobs)
    temp <- B*C
    for(j in 1:m)
      {
        N <- N-A*B*temp
        temp <- rbind(temp[(1+2):(m*2),],temp[1:2,])
        }
    D <- M+N
    
    #repeated measures:
##     if(n!=nobs)D <- apply(cbind(unique(id)),1,function(x,y)
##          {
##            return(apply(cbind(y[,id==x]),1,sum))
##          },y=D)
    if(n!=nobs)D <- t(rowsum(t(D),id))
    
    #normalizing cst:

    X.p <- apply(X,2,function(x,y){return(rep(x,y))},y=m)
    #weights
    if(!is.null(weights))
      {
        X.sqrt.w <- sweep(X,2,sqrt(weights),"*")        
        X.p <- apply(X.sqrt.w,2,function(x,y){return(rep(x,y))},y=m)
      }
    
    cst1 <- apply(phi.p*eps,2,sum)
    cst2 <- apply(phi.p^2,2,sum)
    psi <- apply(eps-phi.p,2,function(x,y){return(rep(x,each=y))},y=2)-matrix(cst1,nrow=m*2,ncol= nobs,byrow=T)
    psi.p <- psi-matrix(cst1-cst2,nrow=m*2,ncol= nobs,byrow=T)
    psi.pp <- apply(eps-phi.p,2,function(x,y){return(rep(x,each=y))},y=2)
    phi.pp <- apply(phi.p,2,function(x,y){return(rep(x,each=y))},y=2)

    alpha <- (X.p*psi)%*%t(X.p*phi.pp)
    block.diag <- rep(list(matrix(1,ncol=2,nrow=2)),m)
    alpha[bdiag(block.diag)==0] <- 0
    
    beta <- (X.p*phi.pp*psi.pp)%*%t(X.p*phi.pp)
    gamma <- (X.p*phi.pp)%*%t(X.p*phi.pp*psi.p)
    pre.normC <- -(alpha-beta-gamma)/n

    normC <- solve(pre.normC)
    #IC
    IC <- normC%*%D
    return(IC)
  }

get.sd <- function(IC,n)
  {
    var.beta <- IC%*%t(IC)/n
    se.beta <- sqrt(diag(var.beta)/n)
    return(matrix(se.beta,ncol=1))
  }

bdiag <- function(x){
     if(!is.list(x)) stop("x not a list")
     n <- length(x)
     if(n==0) return(NULL)
     x <- lapply(x, function(y) if(length(y)) as.matrix(y) else
stop("Zero-length component in x"))
     d <- array(unlist(lapply(x, dim)), c(2, n))
     rr <- d[1,]
     cc <- d[2,]
     rsum <- sum(rr)
     csum <- sum(cc)
     out <- array(0, c(rsum, csum))
     ind <- array(0, c(4, n))
     rcum <- cumsum(rr)
     ccum <- cumsum(cc)
     ind[1,-1] <- rcum[-n]
     ind[2,] <- rcum
     ind[3,-1] <- ccum[-n]
     ind[4,] <- ccum
     imat <- array(1:(rsum * csum), c(rsum, csum))
     iuse <- apply(ind, 2, function(y, imat) imat[(y[1]+1):y[2],
(y[3]+1):y[4]], imat=imat)
     iuse <- as.vector(unlist(iuse))
     out[iuse] <- unlist(x)
     return(out)
}

get.test.stat <- function(coef.n,IC.n,m,n)
  {
    coef.X <- coef.n[extract.X <- seq(2,2*m,by=2),,drop=FALSE]
    
    Sigma.n <- (IC.n%*%t(IC.n)/(n^2))[extract.X,extract.X]
    i.n <- solve(Sigma.n)
    eigen.n <- eigen(i.n)
    T.n <- eigen.n$vectors
    D.n <- diag(eigen.n$values)
    D.sqrt.n <- sqrt(D.n)
    chi.n <- sum((T.n%*%D.sqrt.n%*%t(T.n)%*%coef.X)^2)
    return(chi.n)
  }

get.dsa.fit <- function(dsa.res) {
  model.mats <- get.model.matrices(dsa.res$formula, data = dsa.res$data)
  df <- as.data.frame(cbind(model.mats$Y.orig,model.mats$X))
  wts <- as.numeric(dsa.res$weights[nrow(dsa.res$weights),])
  if (length(wts) == 0)
    wts <- NULL
  model <- dsa.res$model.selected
  environment(model) <- environment()
  if (dsa.res$family == "multinomial")
    fit <- fitModel(model, family = dsa.res$family, data = df, weights = wts)
  else
    fit <- glm(model, family = dsa.res$family, data = df, weights = wts) #avoids creating duplicate methods for glm objects

  return(fit)
}

get.contrast.list <- function(formula, data) {
  mt <- terms(formula, data = data)
  variables <- attr(mt, "variables")
  response <- attr(mt, "response")
  cols <- eval(variables, data)[-response]
  var.names <- as.list(variables)[c(-1,-2)]
  
  lst <- lapply(cols, function(x) {
    if (is.factor(x)) {
      xx <- diag(1, nrow = length(levels(x)))
      colnames(xx) <- rownames(xx) <- as.character(levels(x))
      xx
    }
    else {
      NULL
    }
  })
  names(lst) <- var.names

  ##
  ## remove entries which are not null, no list if
  ## if there are no factors. (there must be a better way)
  ##
  res <- list()
  na <- c()
  if (length(lst) > 0) {
    for (i in 1:length(lst)) {
      if (!is.null(lst[[i]])) {
        res <- c(res, list(lst[[i]]))
        na <- c(na, (as.character(var.names[[i]])))
      }
    }
    names(res) <- na
  }
  if (length(res) <= 0)
    return(NULL)
  else
    return(res)
}

get.model.matrices <- function(formula, data) {
  ## this was lifted from nnet. 
  class.ind <- function(cl) {
    n <- length(cl)
    x <- matrix(0, n, length(levels(cl)))
    x[(1:n) + n * (as.vector(unclass(cl)) - 1)] <- 1
    x[is.na(as.vector(unclass(cl))),] <- NA
    dimnames(x) <- list(names(cl), levels(cl))
    x
  }
  
  complete.formula <- as.formula(paste(deparse(as.list(formula)[[2]]),
                                       " ~ .", sep = ""))
  complete.mf <- model.frame(complete.formula, data = data,
                             na.action = na.pass)
  mt <- attr(complete.mf, "terms")

  ## transpose and drop the intercept column. 
  X <- model.matrix(mt, complete.mf, contrasts.arg =
                    get.contrast.list(complete.formula, data))[,-1,drop=FALSE]

  ## preserve the fact that Ylearn is a matrix, 
  if (is.null(dim(Y <- Y.orig <- model.response(complete.mf)))) {
    Y <- Y.orig <- complete.mf[, 1, drop = FALSE]
  }

  class.labels <- NULL
  if (ncol(Y) == 1) {
    my <- Y[,1]
    if (is.factor(my)) {
      my <- class.ind(my)
      class.labels <- colnames(my)

      ## drop the first column (the class labels will keep the first label):
      my <- my[,-1]
    }
    ## set it back.
    Y <- as.numeric(my)
  }
  ## force it to be a one column matrix.
  Y <- as.matrix(Y)
  colnames(Y.orig) <- colnames(Y.orig) 
  
  return(list(X = X, Y = Y, class.labels = class.labels, Y.orig = Y.orig,
              Y.names = colnames(Y.orig)))
}

get.final.model <- function(modelselected, namesxmat, namesymat)
{
  if(sum(modelselected)==0)res <- paste(namesymat," ~ 1",sep="")
  else {
    res <- apply(modelselected,1,function(x,y) {
      vars <- y[as.logical(x)]
      powers <- x[x!=0]
      termcomponents <- apply(cbind(vars,powers),1,function(x) {
        return(paste(x,collapse="^"))
      })
      term <- ""
      if(length(termcomponents)>0)
        term <- paste(" I(", paste(termcomponents, collapse="*") ,")",sep="")
      return(term)
    }, y=namesxmat)
    if (length(namesymat) == 1) {
      res <- paste(namesymat," ~",paste(res[res!=""],collapse=" +"),sep="")
    }
    else {
      res <- paste("cbind(", namesymat[1], ",", namesymat[2], ")",
                   paste(" ~", paste(res[res!=""],collapse=" +"),sep=""))
    }
  }
  
  ff <- tryCatch(as.formula(res), error = function(...) {
    warning("Unable to parse resulting formula, Potential problem with factor/level conversion.")
    return(res)
  })
  return(ff)
}

##
## formula : A proper formula with factors or with indicators for the factors.
## data :    The data frame
## Xlearn :  The design matrix.
## use.original.df : means that we use the data frame where a factor is one column.
##                   In this case we can get a design matrix using the formula and
##                   the data.frame and we will force in the newly created indicator
##                   columns. This is the case with a model and a data frame directly
##                   from a user. The other case is when we wish to call the DSA using
##                   the original data.frame and a model which the DSA spits back. 
##                   
##
get.forced.terms <- function(formula, data, Xlearn, use.original.df) {
  forcedterms <- NULL

  if (use.original.df) {
    mm <- model.matrix(terms(formula, data = data), data = data)
    
    tryCatch(fitModel(formula,data), error = function(x) {
      stop("Illegal forced terms; Design matrix singular.") #used to try solve(t(mm) %*% mm)
    })
    mm <- mm[, -1, drop=FALSE]
  }
  else {
    df <- get.model.matrices(formula, data)
    df <- as.data.frame(cbind(df$Y.orig, df$X))
    mm <- model.matrix(terms(formula, data = df), data = df)[,-1,drop=FALSE]
  }
  
  Fterms <- colnames(mm)
  nFterms <- length(Fterms)
  
  if (nFterms != 0) {
    forcedterms <- matrix(0, nrow = nFterms, ncol = ncol(Xlearn))
    
    for(i in 1:nFterms) {
      ftermi <- unlist(strsplit(Fterms[i], split="*", fixed=T))
      ftermi <- unlist(apply(cbind(ftermi), 1, function(x)return(strsplit(x,split="I(",fixed=T))))
      ftermi <- unlist(apply(cbind(ftermi), 1, function(x)return(strsplit(x,split=")",fixed=T))))
      ftermi <- unlist(apply(cbind(ftermi), 1, function(x)return(strsplit(x,split=" ",fixed=T))))
      ftermi <- ftermi[ftermi != ""]
      nvari <- length(ftermi)
      powi <- rep(1,nvari)

      for(j in 1:nvari) {
        fvari <- unlist(strsplit(ftermi[j],split="^",fixed=T))
        if(length(fvari)==2) {
          powi[j] <- as.numeric(fvari[2])
          ftermi[j] <- fvari[1]
        }

        if(all(ftermi[j] != colnames(Xlearn)))
          stop("\nWARNING: Incorrect formula for forcedterms.\n")
        
        forcedterms[i, ftermi[j] == colnames(Xlearn)] <- powi[j]
      }
    }
  colnames(forcedterms) <- colnames(Xlearn)
  }

  return(forcedterms)
}


