.libPaths("/work1/Software/R3/lib")
.libPaths()

library(readr)
library(LDheatmap)
library(dplyr)
library(stringr)

library(argparser)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Usage: Rscript script.R --file LD_matrix_file --out output_file")
}

file <- NULL
outfile <- NULL
for (i in seq_along(args)) {
  if (args[i] == "--file" && length(args) > i) {
    file <- args[i + 1]
  } else if (args[i] == "--out" && length(args) > i) {
    outfile <- args[i + 1]
  }
}

if (is.null(file) || is.na(file) || !file.exists(file)) {
  stop("Invalid input file")
} else if (is.null(outfile) || is.na(outfile) || length(outfile) == 0) {
  stop("Invalid output file")
}



LD_file=read_tsv(file)
head(LD_file)[,1:10]
LD=LD_file[,-1]
LD=as.matrix(LD) #必须是matrix，才能进行后续LDheatmap绘图

#获取距离
dis = str_split(colnames(LD),":")
Chr_name=dis[[1]][1]
dis=unlist(lapply(dis,function(x) x[2]))
distance=as.numeric(dis)

#设定颜色
rgb.palette <- colorRampPalette(rev(c("yellow", "red")), space = "rgb")

filename=paste(outfile,sep="")
png(filename,type="cairo",height=1600,width=2400)
LDheatmap(LD, distance,color=rgb.palette(100),flip=TRUE,title=NULL)
dev.off()

pdf(paste0(filename,'.pdf'))
LDheatmap(LD, distance,color=rgb.palette(100),flip=TRUE,title=NULL)
dev.off()

