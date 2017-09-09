#####   Distributed Zero Inflated Poisson Regression model with Gamma Lasso Fit    ######

## define class

setClass("dmrcoef", contains="dgCMatrix")

## A single run implementation of the ZION model

zion_run <- function(xj, argl){
  if(length(xj@i)==0) return(NULL) # n'er occurs
  arglD = arglC = argl
  arglD$family = "binomial"
  arglD$y <- sapply(xj, function(w) if(w > 0) return(1) else return(0))

  if(length(which(arglD$y==0)) == 0){
    fitD <- NULL
    pos <- which(arglD$y==1)
    arglC$family = "poisson"
    arglC$y <- xj[pos]
    arglC$x <- as.matrix(argl$x[pos,])
    if(!is.null(arglC$shift)){arglC$shift <- arglC$shift[pos]}
    if(arglC$cv) {fitC <- do.call(cv.gamlr,arglC)} else {fitC <- do.call(gamlr,arglC)}
  }else if(length(which(arglD$y==1)) == 0){
    fitC <- NULL
    fitD <- NULL
  }else{
    if(arglD$cv) {fitD <- do.call(cv.gamlr,arglD)} else {fitD <- do.call(gamlr,arglD)}
    pos <- which(arglD$y==1)
    arglC$family = "poisson"
    arglC$y <- xj[pos]
    arglC$x <- matrix(argl$x[pos,], nrow=length(pos))
    if(!is.null(argl$shift)){arglC$shift <- argl$shift[pos]}
    if(arglC$cv) {fitC <- do.call(cv.gamlr,arglC)} else {fitC <- do.call(gamlr,arglC)}
  }

  ll <- list("fitD" = fitD, "fitC" = fitC)
  return(ll)
}

## main function - zero inflated Poisson fitting model

zion <- function(cl, counts, covars, mu=NULL, bins=NULL, verb=0, cv=FALSE, ...)
{
  ## check if cl is actually a cluster

  if(!is.null(cl)){
    if(!inherits(cl,"cluster")) stop("first argument `cl' must be NULL or a socket cluster.")
  }


  ## Build an argument list for the gamlr function

  argl <- list()
  argl$family <- NULL
  if(is.null(argl$nlambda))
    argl$nlambda <- formals(gamlr)$nlambda
  argl$verb <- max(verb-1,0)
  argl$cv <- cv
  argl$shift <- mu
  argl$x <- covars

  ## Print the characteristics of the model to be fitted

  cat(sprintf("Fitting the ZION model: %d observations on %d features, %d covariates.\n",
              nrow(counts), ncol(counts), ncol(argl$x)))

  stopifnot( all(counts == floor(counts)) ) ## check if the matrix is an integer matrix

  ## if data frame, how to handle the counts matrix

  if(inherits(counts,"data.frame")){
    if(ncol(counts)>1) counts <- as.matrix(counts)
    else counts <- factor(counts[,1])
  }


  counts=as(counts,"dgCMatrix")  ## convert the counts matrix to sparse matrix
  p <- ncol(counts)
  if(is.null(colnames(counts))) colnames(counts) <- 1:p
  n <- nrow(counts)
  if(n != nrow(argl$x))
    stop("counts and covars have a different number of observations")

  ## The feature names

  C <- ifelse(is.null(cl),Inf,length(cl))
  p <- ncol(counts)
  if(is.null(colnames(counts))) {vars <- 1:dim(counts)[2]} else {vars <- colnames(counts)}
  counts2 <- counts

  ##  Break the counts data into chunks for parallel processing

  if(C < p/4){
    chunks <- round(seq(0,p,length.out=C+1))
    counts <- lapply(1:C,
                     function(i) counts[,(chunks[i]+1):chunks[i+1]])
    counts <- parLapply(cl,
                        counts,
                        function(x)
                          sapply(colnames(x),
                                 function(j) x[,j,drop=FALSE]))
    counts <- unlist(counts,recursive=FALSE)
  } else{
    counts <- sapply(vars,
                     function(j) counts[,j,drop=FALSE]) }


  ##  Parallel implementation of the ZION model - with class specification

  mods <- parLapply(cl,counts, zion_run,argl=argl)

  class(mods) <- "dmr"
  attr(mods,"nobs") <- argl$nobs
  attr(mods,"nlambda") <- argl$nlambda
  attr(mods,"mu") <- argl$shift

  fitC_coef <- as.matrix(do.call(cbind, lapply(mods, function(x) return(coef(x$fitC)))))
  fitD_coef <- as.matrix(do.call(cbind, lapply(mods, function(x) return(coef(x$fitD)))))
  colnames(fitC_coef) =  colnames(fitD_coef) = vars

  ll <- list("model" = mods,
             "fitC_coef" = fitC_coef,
             "fitD_coef" = fitD_coef)

  return(ll)

}

coef.dmr <- function(object, ...){
  B <- lapply(object,coef, ...)
  failures <- sapply(B,is.null)
  if(any(failures)) B[[which(failures)]] <- Matrix(0)
  bx <- unlist(lapply(B,function(b) b@x))
  bi <- unlist(lapply(B,function(b) b@i))
  bp <- c(0,
          cumsum(unlist(lapply(B,function(b) b@p[-1]))))
  Bs <- sparseMatrix(i=bi+1,p=bp,x=bx,
                     dims=c(nrow(B[[1]]),length(B)),
                     dimnames=list(rownames(B[[1]]),names(B)))
  Bs <- as(as(Bs,"dgCMatrix"),"dmrcoef")
  return(Bs)
}
