rm(list = ls())
library(ggplot2)
library(dplyr)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Usage: Rscript script.R --file LD_matrix_file --out output_file --top n")
}

file <- NULL
outfile <- NULL
top_n <- 5
for (i in seq_along(args)) {
  if (args[i] == "--file" && length(args) > i) {
    file <- args[i + 1]
  } else if (args[i] == "--out" && length(args) > i) {
    outfile <- args[i + 1]
  } else if (args[i] == "--top" && length(args) > i) {
    top_n <- as.integer(args[i + 1])
  }
}

if (is.null(file) || is.na(file) || !file.exists(file)) {
  stop("Invalid input file")
} else if (is.null(outfile) || is.na(outfile) || length(outfile) == 0) {
  stop("Invalid output file")
}



#setwd("F:\\01_project\\Bioinfo\\pipline\\GWAS单倍型block分析")
#file = 'match_test.xls'

data = read.table(file,header = T)

# 过滤杂合和缺失的单倍型 提取频数大于20的单倍型的表型
result <- data %>%
  filter(!grepl("H", Geno)) %>%
  filter(!grepl("-", Geno)) %>%
  group_by(Hap) %>%
  mutate(freq = n()) %>%
  ungroup() %>%
  filter(freq > 20) %>%
  select(-freq)

# 判断数组是否为空
if (nrow(result) == 0) {
  # 如果为空，则执行以下代码
  freq_table <- as.data.frame(table(data$Hap))
  freq_table <- freq_table[order(-freq_table$Freq),]
  # 获取前10个最常见的单倍型
  top_types <- head(freq_table$Var1, top_n)
  # 筛选出包含前10个最常见单倍型的行
  result <- data[data$Hap %in% top_types,]
}

p = ggplot(result,aes(x=Hap, y=Pheno,group=Hap)) + 

  geom_boxplot(width = 0.5) +
  #stat_compare_means(comparisons = group_label) +    # 来自R包 ggpubr 用于T检验
  geom_jitter(aes(color=Hap),width = 0.2, alpha = 0.5) + 
  theme_bw()+
  theme(panel.grid=element_blank(),legend.position = "None")
p
ggsave(p,filename=outfile)
