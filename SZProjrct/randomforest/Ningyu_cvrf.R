library(randomForest)
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
trainaccuracyset <- double(nm)
testaccuracyset <- double(nm)

set.seed(1234)
for(j in 1:nm){
  trainaccuracy <- 0
  testaccuracy <- 0
  for (i in 1:10){
    testset <- dataset[c((17 * i -16):(17*i)),]
    trainset <- dataset[-c((17 * i -16):(17*i)),]
    rf <- randomForest(as.factor(state2)~., data = trainset, ntree = 1000, mtry = j)
    trainpredicterror <- predict(rf, trainset, type = 'class')
    testpredicterror <- predict(rf, testset, type = 'class')
    traintable <- table(Predictions=trainpredicterror, actual=trainset$state2)
    testtable <- table(Predictions=testpredicterror, actual=testset$state2)
    trainaccuracy <- trainaccuracy + sum(diag(traintable))/sum(traintable)
    testaccuracy <- testaccuracy + sum(diag(testtable))/sum(testtable)
  }
  trainaccuracyset[j] <- trainaccuracy/10
  testaccuracyset[j] <- testaccuracy/10
}
setwd("C:/Users/Ningyu Sha/Dropbox/Ningyu Sha/Dropbox/Digging-into-metagenomic-data/SZProjrct/randomforest") 
pdf('cross validationn.pdf', width = 10, height = 10)
plot(1:50,testaccuracyset, type = "b", pch = 16, col = c("orange"),xlab = "mtry", ylab = "accuracy" )

#plot(30:500,testaccuracyset, type = "o", xlab = "mtry", ylab = "test set prediction error" )
dev.off()