# _author_ = 'sixueyang'
# _date_ = 2023/5/4 17:00
import pandas as pd
import sys
import os
import numpy as np
import argparse






def change_capitalize(value):
    if isinstance(value,str):
        return value.capitalize()
    else:
        return value

# 负值归0
def trun2zero(value:int):
    if value < 0:
        return 0
    else:
        return value

# 转换整数类型
def trun2int(value):
    return int(value)

# 预处理
def prepare_data(my_model,path,snp_cut:int,flag:int):
    if my_model == 'fst':
        chr_tag = 'CHROM'
        star_tag = 'BIN_START'
        end_tag = 'BIN_END'
        value_tag = 'WEIGHTED_FST'
        SNP_tag = 'N_VARIANTS'

    elif my_model == 'pi':
        chr_tag = 'CHROM'
        star_tag = 'BIN_START'
        end_tag = 'BIN_END'
        value_tag = 'PI'
        SNP_tag = 'N_VARIANTS'

    elif my_model == 'xpclr':
        chr_tag = 'chrom'
        star_tag = 'start'
        end_tag = 'stop'
        value_tag = 'xpclr'
        SNP_tag = 'nSNPs'

    else:
        sys.exit('无该方法，请选择正确的方法！')
    out_value = os.path.basename(path)
    df = pd.read_csv(path,sep='\t').loc[:,[chr_tag,star_tag,end_tag,SNP_tag,value_tag]]
    # print(df)
    df = df[df[SNP_tag] >= snp_cut]
    df[value_tag] = df[value_tag].map(trun2zero)
    df[chr_tag] = df[chr_tag].map(change_capitalize)
    normal_columns = ['CHROM', 'BIN_START', 'BIN_END', f'SNP_Num_{flag}',out_value]
    df.columns = normal_columns
    return df




def setindex(df:pd.DataFrame):
    return df.set_index(['CHROM', 'BIN_START', 'BIN_END'])


# 取并集标签
def merge_dfs(df_list):

    df_list = [i for i in setindex(df_list)]

    merge = pd.concat(df_list,axis=1)
    df_sorted = merge.sort_values(by=['CHROM', 'BIN_START', 'BIN_END'], ascending=True)
    for tag in df_sorted.columns:

        df_sorted[tag] = np.where(df_sorted[tag].isnull(), 0, 1)

    return df_sorted.fillna("NA"),merge.sort_values(by=['CHROM', 'BIN_START', 'BIN_END'], ascending=True).fillna("NA")


# 提取SNP_num相关的列
def get_snp_col(cols,tag):
    snp_lis = []
    for col in cols:
        if tag in col:
            snp_lis.append(col)

    return snp_lis



# 设置参数开关
def get_args():
    parser = argparse.ArgumentParser(description="合并选择分析结果")
    parser.add_argument("-pi1", "--pi1", help="pi结果1", required=True)
    parser.add_argument("-pi2", "--pi2", help="pi结果2", required=True)
    parser.add_argument("-xp1", "--xpclr1", help="xpclr结果1",required=True)
    parser.add_argument("-xp2", "--xpclr2", help="xpclr结果2",required=True)
    parser.add_argument("-fst", "--fst", help="fst结果", required=True)
    parser.add_argument("-t", "--snp_cut", help="窗口内SNP阈值", type=int,default=2)
    parser.add_argument("-o", "--outfile", help="输出文件名", default='my')
    return parser.parse_args()





# 主流程
def main():
    args = get_args()

    fst_file = args.fst
    pi1 = args.pi1
    pi2 = args.pi2
    XPCLR1 = args.xpclr1
    XPCLR2 = args.xpclr2
    snp_cut = int(args.snp_cut)
    out_file = args.outfile


    fst_df = prepare_data('fst',fst_file,snp_cut,1)
    pi1_df = prepare_data('pi',pi1,snp_cut,2)
    pi2_df = prepare_data('pi',pi2,snp_cut,3)
    xp1_df = prepare_data('xpclr',XPCLR1,snp_cut,4)
    xp2_df = prepare_data('xpclr', XPCLR2, snp_cut,5)
    df_list = [fst_df,pi1_df,pi2_df,xp1_df,xp2_df]
    df_list = list(map(setindex,df_list))
    merge_df = pd.concat(df_list,axis=1, join='outer')
    snp_list = get_snp_col(merge_df.columns,'SNP_Num')
    merge_df['SNP_Num'] = merge_df[snp_list].min(axis=1)
    out_col = ['SNP_Num']
    path_list = [fst_file,pi1,pi2,XPCLR1,XPCLR2]
    outval_lis = [os.path.basename(i) for i in path_list]
    out_col.extend(outval_lis)
    out_df = merge_df.loc[:,out_col]
    out_df['SNP_Num'] = out_df['SNP_Num'].map(trun2int)
    # 排序 + 填充NA
    out_sort_df = out_df.sort_index()
    out_sort_df = out_sort_df.fillna('NA')
    out_sort_df.to_csv(out_file, sep='\t', index=True, header=True)



if __name__ == '__main__':
    main()




