#---------------------------------------------------------------------------------------#
# Copyright (c) 2017 Yanmei Ju (BGI-shenzhen). All rights reserved.                      #
# Created by Yanmei Ju (BGI-shenzhen) on 11/11/2017                                     #
# This R program is to calculate the difference of motus in case and control samples    #
# Args:                                                                                 #
#   motu.prof: motus profile, column is sample and row is motu                          #
#   state.prof: state profile, 1st column is sample state: case(1) and control(0)       #
#   prefix: output file prefix                                                          #
# output:                                                                               #
#   5 columns which are estimate; std.error; statistics; pvalue; qvalue respectively    #
# require(pscl) : zero-inflated model
# library(qvalue): adjust p value
#---------------------------------------------------------------------------------------#

# load path
workdir <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(workdir)
rm(list = ls())

require(pscl)
library(qvalue)

# load data
motu.prof <- read.table("../../SZData/mOTU.0.05.profile", header = 1, row.names = 1)
state.prof <- read.table("../../SZData/state171.txt", header = 1, row.names = 1)
prefix <- "zero_infl_wc_state0vs1_2groups_34_1111"
rownames(motu.prof) <- gsub("_", " ", rownames(motu.prof))
colnames(motu.prof) <- gsub("\\.", "-", colnames(motu.prof))

# make motu.prof into integer
motu.prof <- t(motu.prof)
motu.prof <- as.matrix(motu.prof*10^9)

# regress: each motu is y, state (case and control) is x
zero.count <- matrix(0, 360, 1)
regress.result <- matrix(0, 360, 5)
for (i in 1:360) {
  zero.count[i, 1] <- length(which(motu.prof[, i] == 0))
  if(zero.count[i, 1] <= 160) {
    if(zero.count[i, 1] >= 34) {
      zero.infl <- zeroinfl(motu.prof[, i] ~ state.prof$state2)
      zero.infl.summary <- summary(zero.infl)
      regress.result[i, 1:4] <- zero.infl.summary$coefficients$zero[2, ]
    }else if(zero.count[i, 1] < 34) {
      wilc.result <- wilcox.test(motu.prof[, i] ~ state.prof$state2)
      regress.result[i, 1:2] <- c(0, 0)
      regress.result[i, 3] <- wilc.result$statistic
      regress.result[i, 4] <- wilc.result$p.value
    }
  } else {
    regress.result[i, 1:4] <- c(0, 0, 0, 1)
    }
}
rownames(regress.result) <- colnames(motu.prof)

# adjust pvalue of zero.infl and wilcox.test with qvalue method
qvalue.zero.infl.wilcox <- as.data.frame(qvalue(regress.result[which(regress.result[, 3]!=0),4], lambda = seq(0.05, 0.9, 0.05))$qvalues)
qvalue.no.present <- as.data.frame(regress.result[which(regress.result[, 3]==0), 4])
colnames(qvalue.no.present) <- colnames(qvalue.zero.infl.wilcox) <- "qvalue"
qvalue <- rbind(qvalue.zero.infl.wilcox, qvalue.no.present)
regress.qvalue <- qvalue[pmatch(rownames(regress.result), rownames(qvalue)), ]
regress.result[, 5] <- regress.qvalue
colnames(regress.result) <- c("Estimate", "Std.Error", "Zvalue/w", "Pvalue", "Qvalue")

# output: estimate; std.error; statistics; pvalue; qvalue
write.table(regress.result, "zero_infl_wc_state0vs1_2groups_34_1107_bh.txt", quote = F, row.names = T, sep = "\t")




