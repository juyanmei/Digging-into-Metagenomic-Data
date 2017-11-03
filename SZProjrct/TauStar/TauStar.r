
workdir <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(workdir)
rm(list = ls())

#data preprocessing
motu <- read.table("../../SZData/mOTU.0.05.profile",header = 1,row.names = 1)
state <- read.table("../../SZData/state171.txt",header = 1,row.names = 1)
mccb <- read.table("../../SZData/mccb.txt",header = 1,row.names = 1)
mccb <- na.omit(mccb)
panss <- read.table("../../SZData/panss.txt",header = 1,row.names = 1)
mccb <- na.omit(mccb)
colnames(motu) <- gsub("\\.","-",colnames(motu))
mccb_motu <- motu[,pmatch(rownames(mccb),colnames(motu))]
panss_motu <- motu[,pmatch(rownames(panss),colnames(motu))]

calcu_independance <- function(x,y,filename){
pval <- matrix(0,nrow(x),ncol(y))
stat <- matrix(0,nrow(x),ncol(y))
library(TauStar)
for(i in 1:ncol(y)){
  for(j in 1:nrow(x)){
  pval[j,i] = tauStarTest(as.matrix(x[j,]),y[,i])$pVal
  stat[j,i] = tauStarTest(as.matrix(x[j,]),y[,i])$tStar
  rownames(stat) <- rownames(pval) <- rownames(x)
  colnames(stat) <- colnames(pval) <- colnames(y)
  }
}
pvalue <- paste(filename,"pvalue",sep = "_")
write.table(pval,pvalue,quote = F,sep = "\t",row.names = T,col.names = T)
statistic <- paste(filename,"statistic",sep = "_")
write.table(stat,statistic,quote = F,sep = "\t",row.names = T,col.names = T)
}

calcu_independance(mccb_motu,mccb,"mccb_independance")

calcu_independance(panss_motu,panss,"panss_independance")

