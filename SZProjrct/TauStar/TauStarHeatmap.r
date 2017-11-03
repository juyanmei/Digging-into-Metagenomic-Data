workdir <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(workdir)
rm(list = ls())

require(gtools)
data<-read.table("mccb_independance_statistic",row.names=1,header = 1)
sig.data<-read.table("mccb_independance_pvalue",row.names=1,header = 1)

row_id=apply(sig.data,1,function(x){if(any(x<0.05)){a=1}else{a=0}; a})
col_id=apply(sig.data,2,function(x){if(any(x<0.05)){a=1}else{a=0}; a})
data.m <- data[as.logical(row_id),as.logical(col_id)]
sig.data.m <- sig.data[as.logical(row_id),as.logical(col_id)]


library(ComplexHeatmap)
fc_heatmap_t <- Heatmap(
  data.m,
  name = "tStar",
  row_names_side = "right",
  show_row_names = TRUE,
  show_column_names = TRUE,
  row_names_gp = gpar(fontsize = 5),
  column_names_gp = gpar(fontsize = 7),
  # column_names_max_height = unit(5, "cm"),
  # bottom_annotation = ion_bottom,
  # bottom_annotation_height = unit(28, "mm"),
#  width = unit(80, "mm"),
  show_row_dend = TRUE,
  #top_annotation_height = unit(50,"cm"),
  cell_fun = function(j, i, x, y, width, height, fill) {
    if (sig.data.m[i, j] < 0.01 ) {
        grid.text(print("*"),x, y,gp = gpar(fontsize = 4))   #add lable size
    }
    else if (sig.data.m[i, j] < 0.05) {
      grid.text(print("+"), x, y,gp = gpar(fontsize = 4))
    }
  }
  
)
draw(fc_heatmap_t)

