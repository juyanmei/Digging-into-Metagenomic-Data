cv.samQL <- function(x, y, cv.fold, lambda) {
  # cv samQL
  # Args:
  #  x: t(motu) profile
  #  y: phenotype profile
  #  cv.fold: n fold CV
  #  lambda: the lambda wants to used in samQL

  # Returns:
  #  the lambda with min sse
  N <- nrow(x)
  foldid <- sample(rep(seq(cv.fold), length = N))
  samHL.object <- samQL(x, y, lambda = lambda)
  lambda <- samHL.object$lambda
  outlist <- as.list(seq(cv.fold))
  for (i in seq(1, cv.fold)) {
    which <- foldid == i
    if(is.matrix(y)) {
      y.sub <- y[!which, ]
    } else {
      y.sub <- y[!which]
    }
    outlist[[i]] <- samQL(x[!which, ], y.sub, lambda = lambda)
  }
  out <- error.cv(outlist, lambda, x, y, foldid, cv.fold)
  return(out)
}
error.cv <- function(outlist, lambda, x, y, foldid, cv.fold) {
  # the lambda with min sse
  
  # Args:
  #  outlist: the object list of samQL
  #  lambda: samQL's lambda
  #  x: t(motu) profile
  #  y: phenotype profile
  #  foldid: foldid used in cv.samQL
  #  cv.fold: n fold CV used in cv.samQL
  
  # Returns:
  #  the lambda with min sse 
  predmat <- matrix(NA, cv.fold, length(lambda))
  nfolds <- max(foldid)
  for(i in seq(nfolds)) {
    which <- foldid == i
    fitobj <- outlist[[i]]
    preds <- predict(fitobj, x[which, ])$value
    nlami <- length(outlist[[i]]$lambda)
    predmat[i, seq(nlami)] <- apply(preds, 2, function(x) mean((y[which] - x)^2))
  }
  lambda.min <- lambda[which.min(apply(predmat, 2, mean))]
  return(lambda.min)
}