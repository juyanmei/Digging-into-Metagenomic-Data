### canonical correlation analysis
motu <- read.table("../../SZData/mOTU.0.05.profile",header = 1,row.names = 1)
state <- read.table("../../SZData/state171.txt",header = 1,row.names = 1)
patient <- read.table("../../SZData/panss.txt",header = 1,row.names = 1)
motu2 <- as.data.frame(t(motu))

.cran_packages <- c( "shiny","miniUI", "caret", "pls", "e1071", "ggplot2", "randomForest", "dplyr", "ggrepel", "nlme", "devtools",
                     "reshape2", "PMA", "structSSI", "ade4",
                     "ggnetwork", "intergraph", "scales")
.github_packages <- c("jfukuyama/phyloseqGraphTest")
.bioc_packages <- c("genefilter", "impute")
# Install CRAN packages (if not already installed)
.inst <- .cran_packages %in% installed.packages()
if (any(!.inst)){
  install.packages(.cran_packages[!.inst],repos = "http://cran.rstudio.com/")
}
.inst <- .github_packages %in% installed.packages()
if (any(!.inst)){
  devtools::install_github(.github_packages[!.inst])
}

.inst <- .bioc_packages %in% installed.packages()
if(any(!.inst)){
  source("http://bioconductor.org/biocLite.R")
  biocLite(.bioc_packages[!.inst])
}

library(genefilter)
library(PMA)
library(phyloseq)
library(ade4) #dudi.pca
library(ggrepel)

motu2 <- prune_taxa(taxa_sums(motu)>4, motu)
motu2 <- filter_taxa(motu2, filterfun(kOverA(3,2)), TRUE)
X <- otu_table(motu2)
  
ttt = match(rownames(patient),rownames(motu2), nomatch = NA_integer_, incomparables = NULL)
motu2 <- motu2[ttt,]
cca_res2 <- CCA(motu2, patient, penaltyx = .15, penaltyz = .15)
combined <- cbind(motu2[, cca_res2$u != 0],
                  patient[, cca_res2$v != 0])
pca_res <- dudi.pca(combined, scannf = F, nf = 3)
genotype <- substr(rownames(pca_res$li), 1, 2)
sample_type <- substr(rownames(pca_res$l1), 3, 4)
feature_type <- grepl("\\.", colnames(combined))
feature_type <- ifelse(feature_type, "patient", "mOTU")
sample_info <- data.frame(pca_res$li, genotype, sample_type)
feature_info <- data.frame(pca_res$c1,
                           feature = substr(colnames(combined), 1, 6))
ggplot() +  geom_point(data = sample_info,
                       aes(x = Axis1, y = Axis2, col = sample_type, shape = genotype), size = 3) + 
  geom_label_repel(data = feature_info,
                   aes(x = 5.5 * CS1, y = 5.5 * CS2, label = feature, fill = feature_type),
                   size = 2, segment.size = 0.3,
                   label.padding = unit(0.1, "lines"), label.size = 0) +
  geom_point(data = feature_info,
             aes(x = 5.5 * CS1, y = 5.5 * CS2, fill = feature_type),
             size = 1, shape = 23, col = "#383838") +
  scale_color_brewer(palette = "Set2") +
  scale_fill_manual(values = c("#a6d854", "#e78ac3")) +
  guides(fill = guide_legend(override.aes = list(shape = 32, size = 0))) +
  coord_fixed(sqrt(pca_res$eig[2] / pca_res$eig[2])) +
  labs(x = sprintf("Axis1 [%s%% Variance]",
                   100 * round(pca_res$eig[1] / sum(pca_res$eig), 2)),
       y = sprintf("Axis2 [%s%% Variance]",
                   100 * round(pca_res$eig[2] / sum(pca_res$eig), 2)),
       fill = "Feature Type", col = "Sample Type")
# results
# Num non-zeros u's:  10 
# Num non-zeros v's:  1 
# Type of x:  standard 
# Type of z:  standard 
# Penalty for x: L1 bound is  0.15 
# Penalty for z: L1 bound is  0.15 
# Cor(Xu,Zv):  0.6184829
