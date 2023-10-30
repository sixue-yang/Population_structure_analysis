# _author_ = 'sixueyang'
# _date_ = 2023/8/3 17:46
import sys

import pandas as pd
from intervaltree import IntervalTree
from collections import defaultdict
import argparse
import os


# TODO:整合不同流程的top5、注释每个窗口下的基因及功能

# 获取top数据
def get_top(df:pd.DataFrame,chr_tag:str,star_tag:str,end_tag:str,value_tag:str,SNP_tag:str,SNP_cut:int,top_cut:float):
    filter_df = df[df[SNP_tag]>SNP_cut]
    sort_df = filter_df.sort_values(value_tag,ascending=False).head(int(top_cut * filter_df.shape[0]))
    return sort_df.sort_values(by=[chr_tag,star_tag,end_tag])

# 解析bed
def get_tree(bed_path):
    out_data = defaultdict(list)
    with open(bed_path,'r',encoding='utf-8') as file_obj:
        for line in file_obj:
            # tag = line.strip().split('\t')
            chrs,star,end,gene = line.strip().split('\t')
            out_data[chrs].append((int(star),int(end),gene))
    # 构树
    out_data = {key:IntervalTree.from_tuples(out_data[key]) for key in out_data.keys()}
    return out_data

# 解析功能注释
def parse_func(fun_path):
    fun_dic = {}
    with open(fun_path,'r') as file_obj:
        file_obj.readline()
        file_obj.tell()
        for line in file_obj:
            tag = line.strip().split('\t')
            fun_dic[tag[0]] = tag[2]
    return fun_dic

# 注释窗口下的基因
def mk_win_for_gene(df:pd.DataFrame,chr_tag:str,star_tag:str,end_tag:str,bed_tree:dict):
    gene_anno_lis = []
    # func_anno_lis = []
    # chr_group = df.groupby(chr_tag)
    for indexs in df.index:
        chrs = df.loc[indexs,chr_tag]
        star = df.loc[indexs,star_tag]
        end = df.loc[indexs, end_tag]
        if chrs not in bed_tree.keys():
            over_gene_list = ['NA']
        else:
            intervals = bed_tree[chrs].overlap(star,end)
            if len(intervals) == 0:
                over_gene_list = ['NA']
            else:
                over_gene_list = [interval.data for interval in intervals]
        gene_anno_lis.append(';'.join(over_gene_list))
        #func_anno_lis.append(';'.join([fun_dic.get(gene,'') for gene in over_gene_list]))
    anno_df = df.copy()
    anno_df['Gene'] = gene_anno_lis
    # anno_df['Function'] = func_anno_lis
    return anno_df


# 设置参数
def get_args():
    parser = argparse.ArgumentParser(description="选择分析 -- top数据整理.")
    parser.add_argument("-in", "--inputfile", help="选择分析结果", required=True)
    parser.add_argument("-m", "--model", help='分析方法 fst/pi/xpclr',choices=['fst','pi','xpclr'], default='fst')
    parser.add_argument("-b", "--bed", help=r'bed文件路径，文件格式如下，无表头 chrs\tStar\tEnd\tGeneid', required=True)
    #parser.add_argument("-g", "--genename", help=r'基因功能文件 文件格式如下，有表头 Gene_id\tGene_id\tFunction', required=True)
    parser.add_argument("-o", "--out", help='输出文件路径', required=True)
    parser.add_argument("-t", "--topCut", help='选取top的阈值', type=float,default=0.05)
    parser.add_argument("-s", "--SNPCut", help='过滤SNP的阈值', type=int,default=2)

    return parser.parse_args()


# 主函数
def main():
    argv = get_args()
    input_file = argv.inputfile
    my_model = argv.model
    bed_file = argv.bed
    #genename_file = argv.genename
    top = argv.topCut
    snp_cut = argv.SNPCut
    out = argv.out

    in_df = pd.read_csv(input_file,sep='\t',encoding='utf-8')
    bed_tree_dic = get_tree(bed_file)
    #fun_dic = parse_func(genename_file)


    if my_model == 'fst':
        chr_tag = 'CHROM'
        star_tag = 'BIN_START'
        end_tag = 'BIN_END'
        value_tag = 'WEIGHTED_FST'
        SNP_tag = 'N_VARIANTS'
        out_value = 'FST'

    elif my_model == 'pi':
        chr_tag = 'CHROM'
        star_tag = 'BIN_START'
        end_tag = 'BIN_END'
        value_tag = 'PI'
        SNP_tag = 'N_VARIANTS'
        out_value = os.path.basename(input_file).split('.')[0]
    elif my_model == 'xpclr':
        chr_tag = 'chrom'
        star_tag = 'start'
        end_tag = 'stop'
        value_tag = 'xpclr'
        SNP_tag = 'nSNPs'
        out_value = os.path.basename(input_file).split('.')[0]
    else:
        sys.exit('无该方法，请选择正确的方法！')

    com_title = ['CHROM', 'BIN_START', 'BIN_END', 'SNP_Num',out_value]
    # 筛选top5
    top_df = get_top(in_df,chr_tag,star_tag,end_tag,value_tag,SNP_tag,snp_cut,top)
    top_df = top_df.loc[:,[chr_tag,star_tag,end_tag,SNP_tag,value_tag]]
    top_df.columns = com_title

    anno_df = mk_win_for_gene(top_df,'CHROM', 'BIN_START', 'BIN_END',bed_tree_dic)
    # anno_df.columns = com_title
    # anno_df = anno_df.fillna('NA')
    anno_df.to_csv(out,sep='\t',index=False)



if __name__ == '__main__':
    main()



