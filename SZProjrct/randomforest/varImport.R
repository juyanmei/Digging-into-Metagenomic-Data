library(randomForest)
set.seed(1234)
rf12 <- randomForest(as.factor(state2)~., data = dataset, ntree = 500, mtry = 12)

###for Gini
setwd("C:/Users/Ningyu Sha/Dropbox/Ningyu Sha/Dropbox/Digging-into-metagenomic-data/SZProjrct/randomforest") 
pdf('varImpPlot_gini.pdf', width = 10, height = 10)
varImpPlot(rf12, type = 2)
dev.off()
#There are four important bacteria(larger than 1?) how to set threshold

###for OOB
datasetoob = dataset
oob.err = rf12$err.rate
oob.diff = double(360)
for (k in 1:360) {
  permute = rnorm(171, mean=, sd=1)
  datasetoob[, k] = permute
  rfafter <- randomForest(as.factor(state2)~., data = datasetoob, ntree = 500, mtry = 12)
  oob.errafter = rfafter$err.rate
  oob.diff[k] = oob.errafter - oob.err
}