# _author_ = 'sixueyang'
# _date_ = 2023/8/15 15:04
import pandas as pd
import argparse
import os

# TODO:计算两个群体pi之间的比值，群体组合下受选择的pi


def read_pi(pi_path,name):
    df = pd.read_csv(pi_path, sep='\t')
    pi_df = df.loc[:, ['CHROM', 'BIN_START', 'BIN_END','N_VARIANTS', 'PI']]
    pi_df.columns = ['CHROM', 'BIN_START', 'BIN_END',f'N_VARIANTS_{name}', f'PI_{name}']
    return pi_df


# 设置参数
def get_args():
    parser = argparse.ArgumentParser(description="选择分析 -- 不同群体pi 受选择计算.")
    parser.add_argument("-p1", "--pop1", help="群体1 pi值计算结果文件", required=True)
    parser.add_argument("-p2", "--pop2", help="群体2 pi值计算结果文件", required=True)
    parser.add_argument("-o1", "--out1", help='群体1受选择结果', required=True)
    parser.add_argument("-o2", "--out2", help='群体2受选择结果', required=True)
    return parser.parse_args()



def main():
    args = get_args()
    pop1_file = args.pop1
    pop2_file = args.pop2
    out1 = args.out1
    out2 = args.out2

    pop1_name = os.path.basename(pop1_file).split('.')[0]
    pop2_name = os.path.basename(pop2_file).split('.')[0]

    p1_df = read_pi(pop1_file,pop1_name)
    p2_df = read_pi(pop2_file,pop2_name)

    merge_pi = pd.merge(p1_df, p2_df, on=['CHROM', 'BIN_START', 'BIN_END'], how='inner')
    merge_pi['N_VARIANTS'] = [min(i) for i in zip(merge_pi[f'N_VARIANTS_{pop1_name}'].tolist(),merge_pi[f'N_VARIANTS_{pop2_name}'].tolist())]
    # 复制数组，避免后续不同选择结果值发生变化
    p1_select_df = merge_pi.copy()
    p2_select_df = merge_pi.copy()

    p1_select_df['PI'] = merge_pi[f'PI_{pop2_name}']/merge_pi[f'PI_{pop1_name}']
    p2_select_df['PI'] = merge_pi[f'PI_{pop1_name}']/merge_pi[f'PI_{pop2_name}']

    out_p1_select = p1_select_df.loc[:,['CHROM', 'BIN_START', 'BIN_END','N_VARIANTS','PI']]
    out_p2_select = p2_select_df.loc[:, ['CHROM', 'BIN_START', 'BIN_END', 'N_VARIANTS', 'PI']]
    # 输出
    out_p1_select.to_csv(out1,sep='\t',index=False)
    out_p2_select.to_csv(out2, sep='\t', index=False)



if __name__ == '__main__':
    main()




