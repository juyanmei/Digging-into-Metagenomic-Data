motu <- read.table("../../SZData/mOTU.0.05.profile",header = 1,row.names = 1)
state <- read.table("../../SZData/state171.txt",header = 1,row.names = 1)
patient <- read.table("../../SZData/panss.txt",header = 1,row.names = 1)
dpca_re <- ordinate(otu_table(motu, taxa_are_rows, errorIfNULL = TRUE), method = "DPCoA")