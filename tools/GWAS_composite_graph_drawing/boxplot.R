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



#setwd("F:\\01_project\\Bioinfo\\pipline\\GWAS������block����")
#file = 'match_test.xls'

data = read.table(file,header = T)

# �����ӺϺ�ȱʧ�ĵ����� ��ȡƵ������20�ĵ����͵ı���
result <- data %>%
  filter(!grepl("H", Geno)) %>%
  filter(!grepl("-", Geno)) %>%
  group_by(Hap) %>%
  mutate(freq = n()) %>%
  ungroup() %>%
  filter(freq > 20) %>%
  select(-freq)

# �ж������Ƿ�Ϊ��
if (nrow(result) == 0) {
  # ���Ϊ�գ���ִ�����´���
  freq_table <- as.data.frame(table(data$Hap))
  freq_table <- freq_table[order(-freq_table$Freq),]
  # ��ȡǰ10������ĵ�����
  top_types <- head(freq_table$Var1, top_n)
  # ɸѡ������ǰ10����������͵���
  result <- data[data$Hap %in% top_types,]
}

p = ggplot(result,aes(x=Hap, y=Pheno,group=Hap)) + 

  geom_boxplot(width = 0.5) +
  #stat_compare_means(comparisons = group_label) +    # ����R�� ggpubr ����T����
  geom_jitter(aes(color=Hap),width = 0.2, alpha = 0.5) + 
  theme_bw()+
  theme(panel.grid=element_blank(),legend.position = "None")
p
ggsave(p,filename=outfile)