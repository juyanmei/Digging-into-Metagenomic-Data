args <- commandArgs(T)
cat("make sure:\n	y and x must have the same sample id\n")
if(length(args) != 8) {
  stop("Rscript *.R [spe.profile] [env.profile] [step] [nperm] [prefix] [start] [end] [nprocess]")
}
# load package
library(snow)

#
spe.prof <- args[1]
env.prof <- args[2]
penalty.seq <- as.numeric(args[3])
nperms <- as.numeric(args[4])
prefix <- args[5]
start <- as.numeric(args[6])
end <- as.numeric(args[7])
nprocess <- as.numeric(args[8])
# load data
taxo.prof <- read.table(spe.prof, header = 1, row.names = 1)
colnames(taxo.prof) <- gsub("\\.", "-", colnames(taxo.prof))
phe.prof <- read.table(env.prof, header = 1, row.names = 1)
taxo.prof <- taxo.prof[, pmatch(rownames(phe.prof), colnames(taxo.prof))]
taxo.prof <- t(taxo.prof)
phe.prof <- as.matrix(phe.prof)
dim(taxo.prof)
dim(phe.prof)

# tuning para
func <- function(x, taxo.prof, phe.prof, penalty.seq, nperms, start, end) {
  set.seed(0)
  penaltyzs <- x
  penaltyxs <- seq(start, end, penalty.seq)
  cca.perm.out <- matrix(NA, length(penaltyzs), 3)
  for(i in 1:length(penaltyzs)) {
    cca.perm <- CCA.permute(taxo.prof, phe.prof, penaltyxs = penaltyxs, penaltyzs = penaltyzs[i], nperms = nperms)
    cca.perm.out[i, 1] <- cca.perm$zstats[which(penaltyxs == cca.perm$bestpenaltyx)]
    cca.perm.out[i, 2] <- cca.perm$bestpenaltyx
    cca.perm.out[i, 3] <- cca.perm$bestpenaltyz
  }
  cca.perm.out <- as.data.frame(cca.perm.out)
  rownames(cca.perm.out) <- x
  colnames(cca.perm.out) <- c("z", "penaltyx", "penaltyz")
  return(cca.perm.out)
}

mutlinks <- function(cl, taxo.prof, phe.prof, penalty.seq, nperms, start, end) {
  set.seed(0)
  penaltyzs <- seq(0, 1, penalty.seq)
  n <- length(penaltyzs)
  nc <- length(cl)
  options(warn = -1)
  ipenaltyzs <- sapply(split(1:n, 1:nc), function(x) round(x/n, 2))
  options(warn = 0)
  count <- clusterApply(cl, ipenaltyzs, func, taxo.prof, phe.prof, penalty.seq, nperms, start, end)
  count
}

# parallel processing
cl <- makeCluster(nprocess, type="SOCK")
clusterCall(cl, function() { library(PMA); NULL })

#
start
end
#
para.res <- mutlinks(cl, taxo.prof, phe.prof, penalty.seq, nperms, start, end)

penaltyzs <- seq(0, 1, penalty.seq)
res <- matrix(NA, length(penaltyzs), 3)
para.res.nrow.tot <- NULL
for(i in 1:length(para.res)) {
  para.res.nrow <- dim(para.res[[i]])[1]
  para.res.nrow.tot <- cbind(para.res.nrow.tot, para.res.nrow)
}
for(i in 1:length(para.res)) {
  if(i == 1) {
  res[i : (i + para.res.nrow.tot[i] - 1), 1] <- para.res[[i]]$z
  res[i : (i + para.res.nrow.tot[i] - 1), 2] <- para.res[[i]]$penaltyx
  res[i : (i + para.res.nrow.tot[i] - 1), 3] <- para.res[[i]]$penaltyz
  } else {
    res[(sum(para.res.nrow.tot[1:i-1]) + 1) : sum(para.res.nrow.tot[1:i]), 1] <- para.res[[i]]$z
    res[(sum(para.res.nrow.tot[1:i-1]) + 1) : sum(para.res.nrow.tot[1:i]), 2] <- para.res[[i]]$penaltyx
    res[(sum(para.res.nrow.tot[1:i-1]) + 1) : sum(para.res.nrow.tot[1:i]), 3] <- para.res[[i]]$penaltyz 
  }
}

cca.perm.out <- as.data.frame(res)
colnames(cca.perm.out) <- c("z", "penaltyx", "penaltyz")
cca.perm.out$penaltyx <- as.factor(cca.perm.out$penaltyx)
write.table(cca.perm.out, prefix, quote = F, sep = "\t")
#library(ggplot2)
#ggplot(cca.perm.out, aes(penaltyz, z, color = penaltyx)) +
#  geom_point() +
#  facet_grid(. ~ cca.perm.out$penaltyx )






