# _author_ = 'sixueyang'
# _date_ = 2023/8/14 17:30
import pandas as pd
import argparse
import os


# TODO:分割群体文件列表

# 设置参数
def get_args():
    parser = argparse.ArgumentParser(description="选择分析 -- xpclr软件输入vcf文件格式转换.")
    parser.add_argument("-p", "--pop_list", help="原始vcf文件", required=True)
    parser.add_argument("-o", "--outdir", help='输出文件路径', required=True)
    return parser.parse_args()

def main():
    argv = get_args()
    pop_list = argv.pop_list
    out_dir = argv.outdir

    pop_df = pd.read_csv(pop_list,sep='\t',header=None)
    pop_df = pop_df.dropna()
    groups = pop_df.groupby(1)
    for group,group_else in groups:
        out_df = group_else.iloc[:,0]
        out_df.to_csv(out_dir+os.sep+group+'.list',index=False,header=False)


if __name__ == '__main__':
    main()






