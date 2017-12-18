library(randomForest)
set.seed(1234)
p = 171 #number of sample
n = 100 #number of variable
y <- c(1:p)
x <- matrix(0, p, n)
for (ii in 1:p){
  x[ii,] <- rnorm(n, mean = 0, sd = 1)
  if(x[ii,1] <= 0.2){y[ii] <- rbinom(1, 1, 0.2)}
  if(x[ii,1] > 0.2){y[ii] <- rbinom(1, 1, 0.8)}
}

randdataset <- as.data.frame(cbind(x,y))

nm <- 20 #number of mtry
tp <- double(nm)
fp <- double(nm)
oob.err = double(nm)
test.err = double(nm)
trainset <- sample(1:171, 101) #randomly choose 101 train index
for (i in 1:nm){
  set.seed(i)
  rf <- randomForest(as.factor(y)~., data = randdataset, subset = trainset, ntree = 500, mtry = i)#randomForest for classification
  oob.err[i] = rf$err.rate[i,1]
  PredictionWithClass <- predict(rf, randdataset[-trainset, ], type = 'class')
  t <- table(Predictions=PredictionWithClass, actual = (randdataset[-trainset, ]) $ y)
  test.err[i] = sum(diag(t))/sum(t)
  tp[i] <- t[1,1]/(t[1,1]+t[2,1]);
  fp[i] <- t[1,2]/(t[1,2]+t[2,2]);
}

matplot(1:nm, cbind(test.err, oob.err), pch = 19, col = c("red", "blue"), 
        type = "b", ylab = "Error", axes = F)
axis(1)
axis(2)
legend("topright", legend = c("OOB", "Test"), pch = 19, col = c("red", "blue"))

tpp <- tp
fp_list <- sort(fp,index.return= TRUE)
fp <- fp_list$x
ind_fp <- fp_list$ix
temp <- c(1:10)
for (j in 1:10){
  temp[j] <- tp[ind_fp[j]]
}
tp <- temp

#setwd("C:/Users/Ningyu Sha/Dropbox/Ningyu Sha/Dropbox/Digging-into-metagenomic-data/SZProjrct/randomforest") 
#pdf('ROC.pdf', width = 10, height = 10)
#future=paste("ROC curve",".jpg") 
#jpeg(file=future)
#plot(fp,tp, type = "l", xlab = "false positive rate", ylab = "true positive rate" )
#dev.off()