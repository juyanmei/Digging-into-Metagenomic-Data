change <- function(x) {
  # change: change the id from 1 to nrow(x)
  # input:
  #  x: dataframe, y, id, group, time, predictors
  # output:
  #  x: change id form 1 to nrow(x) 
  uniq <- sort(unique(x$id))
  len <- length(uniq)
  for(i in 1:len) {
    x[which(x$id == uniq[i]),"id"] = i
  }
  x
}
occ <- function(x, per) {
  # occ: filter the occurance of bac
  # input:
  #  x: dataframe, y, id, group, time, predictors
  #  per: cutoff of occurance
  # output:
  #  result: after filtering
  
  x.bac <- x[, c(-1:-4)]
  count <- matrix(NA, ncol(x.bac), 1)
  for (i in 1:ncol(x.bac)) {
    count[i, 1] <- length(which(x.bac[, i] != 0))
  }
  rownames(count) <- colnames(x.bac)
  num <- per*nrow(x)*0.01
  x.bac.cut <- x.bac[pmatch(names(count[which(count > num),]), colnames(x.bac))]
  result <- cbind(x[, c(1, 2, 4)], x.bac.cut)
  return(result)
}
outpgee <- function(x) {
  # outpgee: output pgee result
  # input:
  #  x: dataframe, y, id, group, time, predictors
  # output:
  #  list: list[[1]] is cv result; list[[2]] is non-zero variables
  
  test <- x[, c(2,1,3,4:ncol(x))]
  formula <- "y ~.-id"
  family <- gaussian(link = "identity")
  lambda.vec <- seq(0.1,1,0.01)
  set.seed(0)
  cv <- CVfit(formula = formula, id = id, data = test, family = family, scale.fix =
                TRUE,scale.value = 1, fold = 4, lambda.vec = lambda.vec, 
              pindex = c(1,2), eps =10^-6,maxiter = 30, tol = 10^-6)
  myfit1 <- PGEE(formula = formula, id = id, data = test, na.action = NULL,
                 family = family, corstr = "AR-1", Mv = NULL,
                 beta_int = c(rep(0,dim(test)[2]-1)), R = NULL, scale.fix = TRUE,
                 scale.value = 1, lambda = 0.95, pindex = c(1,2), eps = 10^-6,
                 maxiter = 30, tol = 10^-6, silent = TRUE)
  
  index1 <- which(abs(coef(summary(myfit1))[,"Estimate"]) > 10^-3)
  # see the PGEE summary statistics of these non-zero variables
  result <- coef(summary(myfit1))[index1, 1]
  list(cv, result)
}