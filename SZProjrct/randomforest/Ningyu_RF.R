install.packages('randomForest')
install.packages('caTools')
library(randomForest)

# data process, get the training and testing data

motu <- read.table("../../SZData/mOTU.0.05.profile",header = 1,row.names = 1)
state <- read.table("../../SZData/state171.txt",header = 1,row.names = 1)
motu2 <- as.data.frame(t(motu))
dataset <- cbind(motu2,state)
dataset <- dataset[,-c(362)]

names(dataset)[names(dataset)=="Bifidobacterium_catenulatum-Bifidobacterium_pseudocatenulatum_complex"]="BcBp"
names(dataset)[names(dataset)=="Bacteroides_dorei/vulgatus"]="Bdv"
names(dataset)[names(dataset)=="Cronobacter_sakazakii/turicensis"]="Cst"
names(dataset)[names(dataset)=="Bacillus_anthracis/cereus/thuringiensis"]="Bact"
names(dataset)[names(dataset)=="butyrate-producing_bacterium"]="bpb"
names(dataset)[names(dataset)=="[Ruminococcus]_gnavus"]="R_g"
names(dataset)[names(dataset)=="butyrate-producing_bacterium_SS3/4"]="bpb34"
names(dataset)[names(dataset)=="Ruminococcus_sp._SR1/5"]="rsr15"
names(dataset)[names(dataset)=="[Ruminococcus]_obeum"]="R_o"
names(dataset)[names(dataset)=="Clostridium_sp._L2-50"]="Csp50"
names(dataset)[names(dataset)=="[Ruminococcus]_torques"]="R_t"
names(dataset)[names(dataset)=="[Clostridium]_bartlettii"]="C_b"
names(dataset)[names(dataset)=="Enterobacter_hormaechei/cloacae"]="Ehc"
names(dataset)[names(dataset)=="Klebsiella_variicola/pneumoniae"]="Kvp"


nm <- 50 #number of mtry
#tp <- double(nm)
#fp <- double(nm)
oob.err = double(nm)
test.err = double(nm)
trainset <- sample(1:171, 101) #randomly choose 101 train index

for (i in 1:nm){
  set.seed(i)
  rf <- randomForest(as.factor(state2)~., data = dataset, subset = trainset, ntree = 500, mtry = i)
  PredictionWithClass <- predict(rf, dataset[-trainset, ], type = 'class')
  t <- table(Predictions=PredictionWithClass, actual=(dataset[-trainset, ]) $ state2)
  test.err[i] = 1 - sum(diag(t))/sum(t)
  #tp[i] <- t[1,1]/(t[1,1]+t[2,1]);
  #fp[i] <- t[1,2]/(t[1,2]+t[2,2]);
}

setwd("C:/Users/Ningyu Sha/Dropbox/Ningyu Sha/Dropbox/Digging-into-metagenomic-data/SZProjrct/randomforest") 
pdf('error analysis.pdf', width = 10, height = 10)
matplot(1:50, cbind(test.err, oob.err), pch = 19, col = c("red", "blue"), 
        type = "b", ylab = "Error", axes = F)
axis(1)
axis(2)
legend("topright", legend = c("OOB", "Test"), pch = 19, col = c("red", "blue"))
dev.off()

###draw ROC curve
tpp <- tp
fp_list <- sort(fp,index.return= TRUE)
fp <- fp_list$x
ind_fp <- fp_list$ix
temp <- c(1:500)
for (j in 1:500){
  temp[j] <- tp[ind_fp[j]]
}
tp <- temp

setwd("C:/Users/Ningyu Sha/Dropbox/Ningyu Sha/Dropbox/Digging-into-metagenomic-data/SZProjrct/randomforest") 
pdf('ROC.pdf', width = 10, height = 10)
#future=paste("ROC curve",".jpg") 
#jpeg(file=future)
plot(fp,tp, type = "l", xlab = "false positive rate", ylab = "true positive rate" )
dev.off()

#library(pROC)
#auc <- auc(testdata1$state2, PredictionWithProbs[,2])
#plot(roc(testdata1$state2, PredictionWithProbs[,2]))
#importance(rf)
#varImpPlot(rf)