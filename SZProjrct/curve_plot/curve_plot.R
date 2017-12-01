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
# log motu profile
#motu.phe.prof[motu.phe.prof==0] <- 1e-10
#motu.phe.prof <- log10(motu.phe.prof)

ModifymOTU <- function(x, y) {
  # mean motu profile when phe score is same
  #
  # Args:
  #  x: one of the phenotypes of samples
  #  y: motu profile
  #
  # Retures:
  #  mean motu profile when phe score is same
  motu.edit.length <- length(levels(factor(x)))
  count <- 0
  motu.edit.pro <- matrix(0, 360, motu.edit.length)
  rownames(motu.edit.pro) <- rownames(y)
  for (i in 1:motu.edit.length) {
    for (j in 1:360) {
      if (i == 1) {
        motu.edit.pro[j, 1] <- mean(y[j, 1:table(x)[1]])
      } else {
        motu.edit.pro[j, i] <- mean(y[j, (count+1):(count + table(x)[i])])
      }
    }
    count <- count + table(x)[i]
  }
  returnValue(motu.edit.pro)
}

# plot curves im html
for (i in 1:ncol(phe.prof)) {
  phe <- phe.prof[, i]
  names(phe) <- rownames(phe.prof)
  phe <- sort(phe)
  motu <- motu.phe.prof[, pmatch(names(phe), colnames(motu.phe.prof))]
  motu.edit <- ModifymOTU(phe, motu)
  phe.edit <- as.numeric(names(table(phe)))
  iplot.phe <- iplotCurves(motu.edit, phe.edit, 
                           chartOpts = list(xlab = colnames(phe.prof)[i], ylab = "Abundance"))
  iplot.name <- paste(colnames(phe.prof)[i], "html", sep = ".")
  htmlwidgets::saveWidget(iplot.phe, file = iplot.name)
}


for (i in 1:ncol(phe.prof)) {
    phe <- phe.prof[, i]
    names(phe) <- rownames(phe.prof)
    phe <- sort(phe)
    motu <- motu.phe.prof[, pmatch(names(phe), colnames(motu.phe.prof))]
    motu.edit <- ModifymOTU(phe, motu)
    phe.edit <- as.numeric(names(table(phe)))
    iplot.phe <- iplotCurves(motu.edit, phe.edit, 
            chartOpts = list(xlab = colnames(phe.prof)[i], ylab = "Abundance"))
    iplot.name <- paste(colnames(phe.prof)[i], "html", sep = ".")
    htmlwidgets::saveWidget(iplot.phe, file = iplot.name)
}

#heatcols <- heat.colors(100)

pdf('Neg_factor1.pdf', width = 10, height = 4)
heatcols = colorRampPalette(c("red","white","blue"))(100)

j = 1
i = 1
min_val = min(phe.prof[, j])
max_val = max(phe.prof[, j])
plot(x = 1:360, y = motu.phe.prof[ ,1], axes = F, type = 'l', col = heatcols[col_value ], ylim = c(0, 1), 
        ylab = colnames(phe.prof)[1] , xlab = 'motu')
axis(1)
axis(2)

for(i in 2:90){
    col_value = round((phe.prof[i, j] - min_val + 1) / (max_val - min_val) * 100)
    lines(x = 1:360, y = motu.phe.prof[ ,i],  type = 'l', col = heatcols[col_value ])    
}

par(fig=c(.85, .9, .8, .98), mai=c(.1 ,.1 ,.3, .1), new=T)
image(z = 1-rbind(1:length(heatcols)), axes = F, col = heatcols, main = "Weight") #for the legend
box()
axis(2 , at = c(0, 0.25, 0.5, 0.75, 1), labels = c("0", round(sqrt(0.25) * max_val, 1), 
                round(sqrt(0.5) * max_val, 1), round(sqrt(0.75) * max_val, 1) , round(max_val, 1)),las=1)          

dev.off()


