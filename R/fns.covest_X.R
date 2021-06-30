# * Number of Parameters = 1 ----------------------------------------------
#  SINGLE CASE
#' @keywords internal
#' @noRd
thr1.once <- function(X, thr, func){
  output   = list()
  output$S = func(X, thr)
  return(output)
}
#' @keywords internal
#' @noRd
thr1.once.givenS <- function(S, thr, func_S){
  output   = list()
  output$S = func_S(S, thr)
  return(output)
}

#  MULTIPLE CASE
#' @keywords internal
#' @noRd
thr1.singlesum <- function(S1, S2, func_S, thr){
  N = dim(S1)[3]
  output = 0
  for (i in 1:N){
    S1tmp = func_S(S1[,,i], thr)
    S2tmp = S2[,,i]
    output = output + norm(S1tmp-S2tmp,"f")
  }
  return(output)
}


#' @keywords internal
#' @noRd
thr1.multiple_X <- function(X, nCV, nCore, func_X, thrvec, foreach_export = NULL){
  #-------------------------------------------------------------
  # 1. get parameters
  n = nrow(X)
  p = ncol(X)
  # 2. nCV
  n1 = max(round(n*(1-(1/log(n)))), 2)
  # if(n1<2) stop("Sample too small for cross-validation")
  n2 = round(n-n1)
  X_train <-array(0,c(n1,p,nCV))
  S1 = array(0,c(p,p,nCV)) # train
  S2 = array(0,c(p,p,nCV)) # test
  # 3. divide and compute S for nfold cross validation
  idxall = (1:n)
  for (i in 1:nCV){
    idx1 = sample(idxall, n1)
    idx2 = setdiff(idxall, idx1)
    # S1[,,i] = X[idx1,]
    X_train[,,i] = X[idx1,,drop=FALSE]
    S2[,,i] = cov(X[idx2,,drop=FALSE])
  }
  # 4. covariance
  # S = cov(X)

  #-------------------------------------------------------------
  # 1. get ready for multiple thresholding vectors
  stest   = thrvec
  nsearch = length(stest)

  # 2. parallel setup
  cl = parallel::makeCluster(nCore)
  doParallel::registerDoParallel(cl)

  # 3. do parallel computation
  itforeach = NULL
  require(foreach)
  Rs = foreach (itforeach=1:nsearch, .combine=cbind, .export=foreach_export) %dopar% {
    CovTools:::thr1.singlesum(S1=X_train, S2=S2, func_X, stest[itforeach])
  }
  # 4. stop cluster
  parallel::stopCluster(cl)
  # 5. find the minimal element
  thropt = stest[which.min(Rs)]

  #-------------------------------------------------------------
  # 1. compute the optimal one
  outputS = func_X(X, thropt)

  # 2. return output
  output = list()
  output$S  = outputS
  output$CV = data.frame(thr=stest, CVscore=as.vector(Rs))

  # 3. return output
  return(output)
}
