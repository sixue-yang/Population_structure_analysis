# _author_ = 'sixueyang'
# _date_ = 2023/8/16 14:58
import os
import argparse
import pandas as pd

# TODO:合并不同染色体下的xpclr结果
# 对done文件的路径进行解析，获取特定目录下特定受选择方法的所有文件并合并

# 设置参数
def get_args():
    parser = argparse.ArgumentParser(description="选择分析 -- xpclr数据合并.")
    parser.add_argument("-i", "--input", help="初始分析阶段xpclr done文件", required=True)
    parser.add_argument("-o", "--out_prefix", help='输出文件目录+文件头', required=True)
    return parser.parse_args()


# 主程序
def main():
    argv = get_args()
    done_file = argv.input
    out_prefix = argv.out_prefix
    # 提取组合
    data_dir = os.path.dirname(done_file)
    group = os.path.basename(data_dir)

    p1 = group.split('_')[0]
    p2 = group.split('_')[1]

    pop1 = group + '_' + p1
    pop2 = group + '_' + p2
    pop1_files = []
    pop2_files = []
    for file in os.listdir(data_dir):
        if file.startswith(pop1):
            pop1_files.append(os.path.abspath(os.path.join(data_dir,file)))
        elif file.startswith(pop2):
            pop2_files.append(os.path.abspath(os.path.join(data_dir,file)))
        else:
            continue
    pop1_df = pd.concat([pd.read_csv(i,sep='\t',encoding='utf-8') for i in pop1_files],axis=0)
    pop2_df = pd.concat([pd.read_csv(i,sep='\t',encoding='utf-8') for i in pop2_files],axis=0)
    # 输出
    pop1_df.to_csv(out_prefix + f'_{p1}_select.xpclr',sep='\t',index=False)
    pop2_df.to_csv(out_prefix + f'_{p2}_select.xpclr', sep='\t', index=False)



if __name__ == '__main__':
    main()

