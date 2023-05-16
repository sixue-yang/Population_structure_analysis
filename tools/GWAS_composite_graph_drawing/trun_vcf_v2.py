# _author_ = 'sixueyang'
# _date_ = 2023/4/14 15:28
import pandas as pd
import subprocess as sp
import argparse
from multiprocessing import Pool


# 输入vcf 生成ped和map
def mk_ped(vcf_path,out_pre,plink):
    cmd = f'{plink} --vcf {vcf_path} --recode --allow-extra-chr --const-fid --missing-genotype "-" --out {out_pre}'
    sp.call(cmd,shell=True)
    out_ped = out_pre + '.ped'
    out_map = out_pre + '.map'
    return out_ped,out_map

# 工具函数 判断是否一致
def compare_and_remove(char:str):
    if char[0] == char[1]:
        return char[0]
    else:
        return "H"

# 数据预处理
def prepare_data(data):
    sample_data = data[1]
    sample_num = data.shape[0]
    snp_data = data.iloc[:, 6:]
    snp_num = snp_data.shape[1] / 2
    sample_data = sample_data.reset_index(drop=True)
    # 利用reshape将数据变形为2列
    reshaped_df = snp_data.values.reshape(-1, 2)
    # 合并相邻的两列
    merged_col = reshaped_df[:, 0] + reshaped_df[:, 1]
    new_df = pd.DataFrame(merged_col.reshape(int(sample_num), int(snp_num)),
                          columns=[f'snp{i}' for i in range(1, int(snp_num) + 1)])

    # 输出用于绘制LD连锁
    rmfa = pd.DataFrame(snp_data.apply(lambda x: ''.join(x), axis=1),columns=['seq'])
    rmfa = rmfa.reset_index(drop=True)
    rmfa_df = pd.concat([sample_data, rmfa], axis=1)
    rmfa_df.columns = ['sample','seq']


    b_df = new_df.applymap(compare_and_remove)
    # 合并一个样本所有位点
    result_series = pd.DataFrame(b_df.apply(lambda x: ''.join(x), axis=1),columns=['SNP'])

    trun_data = pd.concat([sample_data,result_series],axis=1)



    return trun_data,rmfa_df


def read_ped(ped_path,chunksize:int,max_pool:int,out_file:str):
    with open(f'{out_file}.snp.list','w') as f,open(f'{out_file}.rmfa','w') as rmfa_obj:
        for i, chunk in enumerate(pd.read_csv(ped_path,sep=' ',header=None,low_memory=False,chunksize=chunksize)):
            with Pool(processes=max_pool) as pool:
                result = pool.apply_async(prepare_data, (chunk,))

                chunk,rmfa = result.get()  # 等待子进程结束后进行后续步骤
                print(chunk.to_csv(sep='\t',index=False,header=False),file=f,end='')
                # print(rmfa)
                # print()
                for index in rmfa.index:
                    # print(index)
                    # print(rmfa.loc[index,:].to_list())
                    ids,seq = rmfa.loc[index,:].to_list()

                    print(f'>{ids}\n{seq}',file=rmfa_obj)
                del chunk # 删除chunk以释放内存

                print(f"Processed {i*chunksize} samples.", end='\n')

    print(f"\nFinished processing all samples.")


# 设置参数
def get_args():
    parser = argparse.ArgumentParser(description="转换vcf位点信息")
    parser.add_argument("-v", "--vcf", help="vcf文件路径", required=True)
    parser.add_argument("-o", "--out_prefix", type=str, help="输出样本文件名", required=True)
    parser.add_argument("-t", "--pool", type=int, help="最大进程数，默认同时进行5个进程", default=5)
    parser.add_argument("-c", "--chunksize", type=int, help="一次读取最大条数，默认一次读取10个样本的数据，位点过多时请设置小一些", default=10)
    parser.add_argument("-p", "--plink", type=str, help="plink路径", default='/work1/Users/yangsixue/tools/conda/envs/bio/bin/plink')

    return parser.parse_args()

if __name__ == '__main__':

    # 参数设置
    args = get_args()
    vcf_path = args.vcf
    outname = args.out_prefix
    chunksize = args.chunksize
    max_pool = args.pool
    plink = args.plink
    # 生成ped和map
    out_ped, out_map = mk_ped(vcf_path,outname,plink)
    # out_ped = outname + '.ped'
    read_ped(out_ped,chunksize,max_pool,outname)
