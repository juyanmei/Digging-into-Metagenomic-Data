#-------------------------------------------------------------------------------------#
# Copyright (c) 2017 Yanmei Ju (BGI-shenzhen). Allrights reserved.                    #
# Created by Yanmei Ju (BGI-shenzhen) on 11/12/2017                                   #
# This R program is using to plot curves about motus and PANSS scales.                #
# Args:                                                                               #
#   motu.prof: motus profile, column is sample and row is motu                        #
#   phe.prof: phenotypes profile, column is phenotype and row is sample               # 
#   prefix: output file prefix                                                        #
# output: total PANSS scales curves with matching mOTUs                               # 
# library(qtlcharts)                                                                  #
#-------------------------------------------------------------------------------------#

workdir <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(workdir)
rm(list = ls())

library(qtlcharts)

# load data
motu.prof <- read.table("../../SZData/mOTU.0.05.profile", header = 1, row.names = 1)
phe.prof <- read.table("../../SZData/panss.txt", header = 1, row.names = 1)

# make motus and phenotyoes in the same order
colnames(motu.prof) <- sub("\\.", "-", colnames(motu.prof))
motu.phe.prof <- motu.prof[, pmatch(rownames(phe.prof), colnames(motu.prof))]
motu.phe.prof <- as.matrix(motu.phe.prof)
#motu.phe.prof <- log10(motu.phe.prof)

# plot curves im html
for (i in 1:ncol(phe.prof)) {
  phe <- phe.prof[, i]
  names(phe) <- rownames(phe.prof)
  phe <- sort(phe)
  motu <- motu.phe.prof[, pmatch(names(phe), colnames(motu.phe.prof))]
  iplot.phe <- iplotCurves(motu, phe, chartOpts = list(xlab = colnames(phe.prof)[i],
                                                       ylab = "Abundance"))
  iplot.name <- paste(colnames(phe.prof)[i], "html", sep = ".")
  htmlwidgets::saveWidget(iplot.phe, file = iplot.name)
}



