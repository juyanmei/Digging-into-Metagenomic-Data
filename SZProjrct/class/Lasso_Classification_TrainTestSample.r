                                                 workdir <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(workdir)
rm(list = ls())

#data preprocessing

motu <- read.table("../../SZData/mOTU.0.05.profile",header = 1,row.names = 1)
state <- read.table("../../SZData/state171.txt",header = 1,row.names = 1)
colnames(motu) <- gsub("\\.","-",colnames(motu))
control_motu <- motu[,which(state$state2==0)]
set.seed(0)
train <- sample(1:81,44)
train_control <- control_motu[,train]
test <- setdiff(1:81,train)
test_control <- control_motu[,test]
train_case <- motu[,which(state$state3==1)]
test_case <- motu[,which(state$state3==2)]
train.x <- cbind(train_control,train_case)
test.x <- cbind(test_control,test_case)
train.y <- state[pmatch(colnames(train.x),rownames(state)),]
test.y <-  state[pmatch(colnames(test.x),rownames(state)),]

#classification

library(glmnet)
source("ROC.r")
grid=10^ seq (10,-2, length =1000)
set.seed(17)
cv.out <- cv.glmnet(t(train.x),train.y[,1],family = "binomial", type.measure="class", nfolds = 5)
lasso.pred <- predict(cv.out$glmnet.fit,s=cv.out$lambda.1se,newx = t(test.x),type="class")
table(lasso.pred[,1],test.y[,1])
plot_roc(lasso.pred[,1],test.y[,1])
lasso.coef <- predict(cv.out$glmnet.fit,type='coefficients',s=cv.out$lambda.1se)[1:ncol(t(test.x)),]
lasso.f <- lasso.coef[lasso.coef!=0]
write.table(lasso.f,"Lasso_Classification_TrainTestSample_result",quote=F,row.names = T,sep = "\t")