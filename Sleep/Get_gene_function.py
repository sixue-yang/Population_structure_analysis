# _author_ = 'sixueyang'
# _date_ = 2023/8/14 15:36
import pandas as pd
import argparse

# TODO:获取方法经过top筛选后的基因及功能

# 解析基因功能
def parse_func(fun_path):
    fun_dic = {}
    with open(fun_path,'r') as file_obj:
        file_obj.readline()
        file_obj.tell()
        for line in file_obj:
            tag = line.strip().split('\t')
            fun_dic[tag[0]] = tag[2]
    return fun_dic



# 设置参数
def get_args():
    parser = argparse.ArgumentParser(description="选择分析 -- top 基因及功能整理.")
    parser.add_argument("-t", "--topfile", help="筛选后的文件", required=True)
    parser.add_argument("-g", "--genename", help=r'基因功能文件 文件格式如下，有表头 Gene_id\tGene_id\tFunction', required=True)
    parser.add_argument("-o", "--out", help='输出文件路径', required=True)
    return parser.parse_args()

def main():
    argv = get_args()
    top_file = argv.topfile
    genename_file = argv.genename
    out_file = argv.out
    genename_dic = parse_func(genename_file)
    out_df = pd.DataFrame()

    top_df = pd.read_csv(top_file,sep='\t',encoding='utf-8')
    top_df = top_df.fillna('NA')
    gene_list = top_df['Gene'].tolist()
    gene_set = set()
    for genes in gene_list:
        # print(genes)
        if genes == 'NA':
            continue
        for gene in genes.split(';'):
            gene_set.add(gene)

    out_df['Gene'] = list(gene_set)
    out_df['Function'] = [genename_dic.get(gene,'NA') for gene in gene_set]

    out_df.sort_values(by='Gene')
    out_df.to_csv(out_file,index=False,sep='\t')



if __name__ == '__main__':
    main()

