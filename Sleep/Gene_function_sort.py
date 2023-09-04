# _author_ = 'sixueyang'
# _date_ = 2023/8/25 9:52
import pandas as pd
import os
import numpy as np
import argparse
# TODO:整理各个方法筛选后基因及对应功能注释文件

# 解析功能文件
def parse_func(function_path):
    fun_dic = {}
    with open(function_path) as fun_obj:
        for line in fun_obj:
            gene_id,func = line.strip().split('\t')
            fun_dic[gene_id] = func
    return fun_dic

def mk_func_dic(gene_path_list:list):
    fuc_dic_list = [parse_func(i) for i in gene_path_list]
    fun_dic = fuc_dic_list[0]
    for i in fuc_dic_list[1:]:
        fun_dic.update(i)
    return fun_dic

def setindex(df_list):

    for df in df_list:
        tmp = df.set_index(['Gene'])
        yield tmp

def prepare(df_path,name):

    df = pd.read_csv(df_path,sep='\t')
    my_df = pd.DataFrame()
    my_df['Gene'] = df.loc[:,'Gene']
    my_df[name] = [1 for i in df.index]

    return my_df


def merge_dfs(df_list):

    df_list = [i for i in setindex(df_list)]

    merge = pd.concat(df_list,axis=1)
    df_sorted = merge.sort_values(by=['Gene'], ascending=True)
    for tag in df_sorted.columns:
        df_sorted[tag] = np.where(df_sorted[tag].isnull(), 0, 1)
    return df_sorted


def pre_name(path):
    return os.path.basename(path)


# 设置参数
def get_args():
    parser = argparse.ArgumentParser(description="选择分析 -- 结果整理 整合一个组合不同方法下筛选的基因.")
    parser.add_argument("-f", "--fst_top", help="fst 筛选的基因文件", required=True)
    parser.add_argument("-p1", "--pi_top1", help='pi 1 筛选的基因文件', required=True)
    parser.add_argument("-p2", "--pi_top2", help='pi 2 筛选的基因文件', required=True)
    parser.add_argument("-x1", "--xpclr_top1", help='XPCLR 1 筛选的基因文件', required=True)
    parser.add_argument("-x2", "--xpclr_top2", help='XPCLR 2 筛选的基因文件', required=True)
    parser.add_argument("-o", "--out_sort", help='整理后文件输出路径', required=True)
    return parser.parse_args()


if __name__ == '__main__':

    # 设置参数
    argv = get_args()
    fst_file = argv.fst_top
    p1_file = argv.pi_top1
    p2_file = argv.pi_top2
    xpclr1_file = argv.xpclr_top1
    xpclr2_file = argv.xpclr_top2
    out_file = argv.out_sort



    path_list = [fst_file,p1_file,p2_file,xpclr1_file,xpclr2_file]
    func_dic = mk_func_dic(path_list)
    name_list = list(map(pre_name,path_list))
    df_list = [prepare(i, j) for i, j in zip(path_list, name_list)]
    all_df = merge_dfs(df_list)
    all_df['Function'] = all_df.index.map(func_dic)
    all_df.fillna('NA')
    all_df.to_csv(out_file, sep='\t')





