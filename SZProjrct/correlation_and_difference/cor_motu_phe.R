#-------------------------------------------------------------------------------------#
# Copyright (c) 2017 Yanmei Ju (BGI-shenzhen). Allrights reserved.                    #
# Created by Yanmei Ju (BGI-shenzhen) on 11/10/2017                                   #
# This R program is using to calculate the correlation between motus and phenotypes.  #
# Args:                                                                               #
#   motu.prof: motus profile, column is sample and row is motu                        #
#   phe.prof: phenotypes profile, column is phenotype and row is sample               #　　　　　
#   cor.method: correlation method, kendall or spearman                               #
#   value.rate: the value rate of motus in samples(%)                                 #
#   prefix: output file prefix                                                        #
# output: correlation and pvalue between motus and phenotypes                         # 
# Don't need third party R package                                                    #
#-------------------------------------------------------------------------------------#

# load path
workdir <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(workdir)
rm(list = ls())

# load data
motu.prof <- read.table("../../SZData/mOTU.0.05.profile", header = 1, row.names = 1)
phe.prof <- read.table("../../SZData/panss.txt", header = 1, row.names = 1)
cor.method <- "kendall"
value.rate <- 20
prefix <- "cor_1110"

# make motus and phenotyoes in the same order
colnames(motu.prof) <- sub("\\.", "-", colnames(motu.prof))
motu.phe.prof <- motu.prof[, pmatch(rownames(phe.prof), colnames(motu.prof))]
motu.phe.prof <- as.matrix(motu.phe.prof)

# calculate correlation with specific value rate
estimate <- matrix(0, nrow(motu.phe.prof), ncol(phe.prof))
pvalue <- matrix(0, nrow(motu.phe.prof), ncol(phe.prof))
zero.count <- matrix(0,360,1)
value.count <- (1 - value.rate * 0.01) * nrow(phe.prof)
for (i in 1:360) {
  zero.count[i,1] <- length(which(motu.phe.prof[i, ] == 0 ))
}
for (i in 1:nrow(motu.phe.prof)) {
  if (zero.count[i,] <= value.count) {
    for (j in 1:ncol(phe.prof)) {
      estimate[i,j] <- cor.test(motu.phe.prof[i, ], phe.prof[, j], method = cor.method)$estimate
      pvalue[i,j] <- cor.test(motu.phe.prof[i, ], phe.prof[, j], method = cor.method)$p.value 
    }
  } else {
    estimate[i, ] <- rep(0,ncol(phe.prof))
    pvalue[i, ] <- rep(1,ncol(phe.prof))
  }
}
colnames(estimate) <- colnames(pvalue) <- colnames(phe.prof)
rownames(estimate) <- rownames(pvalue) <- gsub("_", " ", rownames(motu.phe.prof))

# output estimate and p.value
filename.estimate <- paste(prefix, "estimate", value.rate, cor.method, sep = "_")
filename.pvalue <- paste(prefix, "pvalue", value.rate, cor.method, sep = "_")
write.table(estimate, file=filename.estimate, append = F, quote = F, sep = "\t", row.names = T)
write.table(pvalue, file= filename.pvalue, append = F, quote = F, sep = "\t", row.names = T)
