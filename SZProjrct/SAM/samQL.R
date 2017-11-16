#-------------------------------------------------------------------------------------#
# Copyright (c) 2017 Yanmei Ju (BGI-shenzhen). Allrights reserved.                    #
# Created by Yanmei Ju (BGI-shenzhen) on 11/16/2017                                   #
# This R program is using to regression                                               #
# Args:                                                                               #
#   motu.prof: motus profile, column is sample and row is motu                        #
#   phe.prof: phenotypes profile, column is phenotype and row is sample               # 
#   prefix: output file prefix                                                        #
# output: total PANSS scales curves with matching mOTUs                               # 
# library(qtlcharts)                                                                  #
#-------------------------------------------------------------------------------------#

#load path
workdir <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(workdir)
rm(list = ls())

# data preprocessing
motu.prof <- read.table("../../SZData/mOTU.0.05.profile", header = 1, row.names = 1)
phe.prof <- read.table("../../SZData/panss.txt", header = 1, row.names = 1)
phe.prof <- na.omit(phe.prof)
colnames(motu.prof) <- gsub("\\.", "-", colnames(motu.prof))
motu.phe <- motu.prof[, pmatch(rownames(phe.prof), colnames(motu.prof))]

# filter occurrence rate
zero.count <- matrix(0, nrow(motu.phe), 1)
rownames(zero.count) <- rownames(motu.phe)
for (i in 1:nrow(motu.phe)) {
  zero.count[i, 1] <- length(which(motu.phe[i, ] == 0))
}
zero.occu <- zero.count[which(zero.count < 73), 1]
motu.reg <- motu.phe[pmatch(names(zero.occu), rownames(motu.phe)), ]
motu.reg <- t(motu.reg)

# train : test set = 7 : 3
set.seed(0)
train <- sample(1:90, 63)
motu.reg.train <- motu.reg[train, ]
phe.prof.train <- phe.prof[train, ]
test <- setdiff(1:90, train)
motu.reg.test <- motu.reg[test, ]
phe.prof.test <- phe.prof[test, ]

# train set : total motu and total phenotype
# test set : total motu and total phenotype
# motu.reg.test <- motu.reg.train <- motu.reg
# phe.prof.test <- phe.prof.train <- phe.prof

# regress
library(SAM)
source("cv.samQL.R")
grid <- 10^ seq (10, -2, length =100)
reg.predict <- matrix(NA, nrow(phe.prof.test), ncol(phe.prof))
for (i in 1:ncol(phe.prof)) {
  set.seed(1)
  phe.reg.train <- phe.prof.train[, i]
  phe.reg.test <- phe.prof.test[, i]
  lambda.min <- cv.samQL(motu.reg.train, phe.reg.train, cv.fold = 5, lambda = grid)
  reg.lambda.min <- samQL(motu.reg.train, phe.reg.train, lambda = lambda.min)
  reg.predict[, i] <- predict(reg.lambda.min, motu.reg.test)$values
}

