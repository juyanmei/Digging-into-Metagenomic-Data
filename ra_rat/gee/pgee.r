#-------------------------------------------------------------------------------#
# Copyright (c) 2018 Yanmei Ju (BGI-shenzhen). Allrights reserved.              #
# Created by Yanmei Ju (BGI-shenzhen) on 01/15/2018                             #
# This R program is using to do sparse gee                                      #
# Args:                                                                         #
#   motu.prof: column is sample, row is motu                                    #
#   per: occurance cutoff 
# output:                                                                       #
#   out: list such as cv, select non-zero variables                             #   
# library(PGEE)                                                                 #
#-------------------------------------------------------------------------------#

# load dir
workdir <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(workdir)
rm(list = ls())

# load packages
library(PGEE)
source("pgee_func.R")

# load data
motu.prof <- read.table("metaphlan_245_rat_total_info.txt", header = 1, row.names = 1)
per <- 0

# modify data
colnames(motu.prof)[1:4] <- c("y", "id", "group", "time") 
gjk <- motu.prof[which(motu.prof$group == 'GJK'), ]
hfd <- motu.prof[which(motu.prof$group == 'HFD'), ]
lcasei <- motu.prof[which(motu.prof$group == 'Lcasei'), ]
model <- motu.prof[which(motu.prof$group == 'Model'), ]
mtx <- motu.prof[which(motu.prof$group == 'MTX'), ]
zqftn <- motu.prof[which(motu.prof$group == 'ZQFTN'), ]

# combine group
model.gjk <- rbind(model, gjk)
model.lcasei <- rbind(model, lcasei)
model.mtx <- rbind(model, mtx)

# order id
model.mtx <- change(model.mtx)
model.lcasei <- change(model.lcasei)
model.gjk <- change(model.gjk)

# filter occurance
model.mtx.occ <- occ(model.mtx, per)
model.lcasei.occ <- occ(model.lcasei, per)
model.gjk.occ <- occ(model.gjk, per)

# pgee
model.mtx.occ.res <- outpgee(model.mtx.occ)
model.lcasei.occ.res <- outpgee(model.lcasei.occ)
model.gjk.occ.res <- outpgee(model.gjk.occ)

# output result
write.table(model.mtx.occ.res[[2]], "model.mtx.pgee.AR1.txt", quote = F, sep = "\t")
write.table(model.lcasei.occ.res[[2]], "model.lcasei.pgee.AR1.txt", quote = F, sep = "\t")
write.table(model.gjk.occ.res[[2]], "model.gjk.pgee.AR1.txt", quote = F, sep = "\t")
