# _author_ = 'sixueyang'
# _date_ = 2023/8/31 13:46
import pandas as pd
import sys
import os
import numpy as np
import argparse
import copy



# 设置传参
def get_args():
    parser = argparse.ArgumentParser(description="选择分析 -- 结果整理 整合一个组合不同方法下筛选的结果.")
    parser.add_argument("-f", "--fst_top", help="fst top文件", required=True)
    parser.add_argument("-p1", "--pi_top1", help='pi 1 top文件', required=True)
    parser.add_argument("-p2", "--pi_top2", help='pi 2 top文件', required=True)
    parser.add_argument("-xp1", "--xpclr_top1", help='XPCLR 1 top文件', required=True)
    parser.add_argument("-xp2", "--xpclr_top2", help='XPCLR 2 top文件', required=True)
    parser.add_argument("-o1", "--out_value", help='窗口对应基因值文件', required=True)
    parser.add_argument("-o2", "--out_sleep_stat", help='窗口受到选择文件', required=True)
    return parser.parse_args()


def setindex(df:pd.DataFrame):
    return df.set_index(['CHROM', 'BIN_START', 'BIN_END'])

def change_capitalize(value):
    if isinstance(value,str):
        return value.capitalize()
    else:
        return value

# 数据框前处理
def read_top(data_path):
    df = pd.read_csv(data_path,sep='\t')
    raw_col = df.columns
    new_col = ['CHROM', 'BIN_START', 'BIN_END']
    new_col.append(raw_col[4])
    read_df = df.loc[:,new_col]
    new_value_name = os.path.basename(data_path)
    _new_col = copy.deepcopy(new_col)
    _new_col[-1] = new_value_name
    read_df.columns = _new_col
    read_df['CHROM'] = read_df['CHROM'].map(change_capitalize)
    return setindex(read_df)

# 合并
def merge_dfs(df_list):


    merge = pd.concat(df_list,axis=1)
    df_sorted = merge.sort_values(by=['CHROM', 'BIN_START', 'BIN_END'], ascending=True)
    for tag in df_sorted.columns:
        df_sorted[tag] = np.where(df_sorted[tag].isnull(), 'No', 'Yes')

    return df_sorted.fillna("NA"),merge.sort_values(by=['CHROM', 'BIN_START', 'BIN_END'], ascending=True).fillna("NA")


# 检查路径
def check_dir(path):
    dir_name = os.path.dirname(path)
    os.makedirs(dir_name, exist_ok=True)


if __name__ == '__main__':
    argv = get_args()
    fst_file = argv.fst_top
    p1_file = argv.pi_top1
    p2_file = argv.pi_top2
    xpclr1_file = argv.xpclr_top1
    xpclr2_file = argv.xpclr_top2
    out_file1 = argv.out_value
    out_file2 = argv.out_sleep_stat
    # 检查目录
    check_dir(out_file1)
    check_dir(out_file2)

    path_lis = [fst_file,p1_file,p2_file,xpclr1_file,xpclr1_file,xpclr2_file]
    df_list = [read_top(i) for i in path_lis]
    win_sleep_stat_df,win_sleep_value_df = merge_dfs(df_list)
    win_sleep_value_df.to_csv(out_file1,sep='\t',index=True,header=True)
    win_sleep_stat_df.to_csv(out_file2,sep='\t',index=True,header=True)



