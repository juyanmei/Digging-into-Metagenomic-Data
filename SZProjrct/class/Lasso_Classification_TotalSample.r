workdir <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(workdir)
rm(list = ls())

#data preprocessing
train.x <- read.table("../../SZData/mOTU.0.05.profile",header = 1,row.names = 1)
train.y <- read.table("../../SZData/state171.txt",header = 1,row.names = 1)
colnames(train.x) <- gsub("\\.","-",colnames(train.x))
train.y <- train.y[pmatch(colnames(train.x),rownames(train.y)),]

#classification
library(glmnet)
grid=10^ seq (10,-2, length =1000)
set.seed(0)
cv.out <- cv.glmnet(t(train.x),train.y[,1],family = "binomial", type.measure="class", nfolds = 5,lambda = grid)
lasso.pred <- predict(cv.out$glmnet.fit,s=cv.out$lambda.1se,newx = t(train.x),type="class")
table(lasso.pred,train.y[,1])
source("ROC.r")
roc1 <- plot_roc(lasso.pred[,1],train.y[,1])
lasso.coef <- predict(cv.out$glmnet.fit,type='coefficients',s=cv.out$lambda.1se)[1:ncol(t(train.x)),]
lasso.f <- lasso.coef[lasso.coef!=0]
write.table(lasso.f,"Lasso_Classification_TotalSample_result",quote=F,row.names = T,sep = "\t")




