
workdir <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(workdir)
rm(list = ls())

#data preprocessing
motu <- read.table("../../SZData/mOTU.0.05.profile",header = 1,row.names = 1)
state <- read.table("../../SZData/state171.txt",header = 1,row.names = 1)
phe <- read.table("../../SZData/panss.txt",header = 1,row.names = 1)
phe <- na.omit(phe)
colnames(motu) <- gsub("\\.","-",colnames(motu))
motu <- motu[,pmatch(rownames(phe),colnames(motu))]
state <- state[pmatch(rownames(phe),rownames(state)),]

train.x <- t(motu)
train.y <- phe
test.x  <- t(motu)
test.y <- phe

#lasso regression, alpha = 1

library(glmnet)
set.seed(0)
grid=10^ seq (10,-2, length =1000)
pseR2 <- matrix(0,ncol(phe),1)
lasso.pred <- matrix(0,nrow(test.y),ncol(phe))
for(i in 1:ncol(phe)){
  set.seed(1)
  cv.out <- cv.glmnet(train.x,train.y[,i],family = "gaussian", alpha = 0.15,lambda = grid,nfolds = nrow(motu))
  bestlam <- cv.out$lambda.min
  lasso.pred[,i] <- predict(cv.out$glmnet.fit,s=bestlam,newx = test.x,type="response")
  pseR2[i,1] <- 1-mean((lasso.pred[,i]-test.y[,i])^2)/mean((test.y[,i]-mean(test.y[,i]))^2)
  rownames(pseR2) <- colnames(phe)
  if(pseR2[i,1]>0.15){
    lasso.coef <- predict(cv.out$glmnet.fit,type='coefficients',s=bestlam)[1:ncol(test.x),]
    lasso.f <- lasso.coef[lasso.coef!=0]
    name <- paste(colnames(test.y)[i],"panss_lasso.txt",sep = "_")
    write.table(lasso.f,name,quote=F,row.names = T,sep = "\t")
  }
}



