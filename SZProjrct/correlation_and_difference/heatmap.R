#---------------------------------------------------------------------------------------#
# Copyright (c) 2017 Yanmei Ju & Jie Zhu(BGI-shenzhen). Allrights reserved.             #
# Created by Yanmei Ju (BGI-shenzhen) on 11/10/2017                                     #
# This R program is using to draw heatmap of significant correlation of enriched motus  #
# Args:                                                                                 #
#   estimate: estimation of correlation                                                 #
#   pvalue: pvalue of estimate                                                          #
#   enrich: significant difference of case and control                                  #
#   prefix: output file prefix                                                          #
# output:                                                                               #
#   heatmap of correlation: pink letters represent motu enrich in case, blue in control # 
# libary(readr): need to load tidyr library for reading files                           #
# library(ComplexHeatmap): to draw heatmap                                              #
# require(circlize)                                                                     #
#---------------------------------------------------------------------------------------#


# load path
workdir <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(workdir)
rm(list = ls())

library(readr)

# load data
estimate <- read_delim("cor_1110_estimate_20_kendall", delim = "\t", col_names = TRUE)
pvalue <- read_delim("cor_1110_pvalue_20_kendall", delim = "\t", col_names = TRUE)
enrich <- read_delim("../../SZData/motu.0.05.enrich", delim="\t", col_names = TRUE)
prefix <- "cor_1110_20_kendall"

# preprocess data
estimate <- estimate[pmatch(as.matrix(enrich[, 1]), as.matrix(estimate[, 1])), ]
pvalue <- pvalue[pmatch(as.matrix(enrich[, 1]), as.matrix(pvalue[, 1])), ]
estimate[is.na(estimate)] <- 0
pvalue[is.na(pvalue)] <- 1

# just draw significant correlation
row.id <- apply(pvalue,1,function(x) {if (any(x < 0.05)) {a = 1} else {a = 0}; a})
col.id <- apply(pvalue,2,function(x) {if (any(x < 0.05)) {a = 1} else {a = 0}; a})
col.id[1] <- 1
estimate.sig <- estimate[as.logical(row.id), as.logical(col.id)]
pvalue.sig <- pvalue[as.logical(row.id), as.logical(col.id)]
enrich.sig <- enrich[pmatch(as.matrix(estimate.sig[,1]), as.matrix(enrich[, 1])), ]

# letter color
enrich.col <- as.matrix(enrich.sig[, 4])
enrich.col[enrich.col == 0] <- "#54AAE9"
enrich.col[enrich.col == 1] <- "#E74899"

# draw data tidy
estimate.draw <- as.data.frame(estimate.sig)
rownames(estimate.draw) <- estimate.draw[, 1]
estimate.draw <- estimate.draw[, -1]
pvalue.draw <- as.data.frame(pvalue.sig)
rownames(pvalue.draw) <- pvalue.draw[, 1]
pvalue.draw <- pvalue.draw[, -1]

# draw
require(circlize)
library(ComplexHeatmap)
fc_heatmap_t <- Heatmap(
  estimate.draw,
  name = "correlation",
  row_names_side = "right",
  row_names_gp = gpar(fontsize = 7,col = enrich.col),
  col = colorRamp2(c(-0.4, 0, 0.4), c("#65A5DA", "white", "#EB3C96")),
  show_row_names = TRUE,
  show_column_names = TRUE,
  column_names_gp = gpar(fontsize = 7),
  # width = unit(50, "mm"),
  show_row_dend = TRUE,
  cell_fun = function(j, i, x, y, width, height, fill) {
    # lable the significant degree, 
    # "*" represents pvalue < 0.01, "+" represents pvalue < 0.05
    
    # Args:
    #   j: column number of pvalue.draw
    #   i: row number of pvalue.draw
    #   x, y, width, height, fill: just use as default
    
    if (pvalue.draw[i, j] < 0.01 ) {
      grid.text(print("*"), x, y, gp = gpar(fontsize = 8))
    }
    else if (pvalue.draw[i, j] < 0.05) {
      grid.text(print("+"), x, y, gp = gpar(fontsize = 8))}
  }
)
fig.name <- paste(prefix, "pdf", sep = ".")
pdf(fig.name)
draw(fc_heatmap_t)
dev.off()

