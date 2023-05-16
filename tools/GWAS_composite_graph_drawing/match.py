# _author_ = 'sixueyang'
# _date_ = 2023/4/17 14:25
import argparse
import sys
import pandas as pd


# 解析单倍型文件 {基因型:单倍型名称}
def read_nex(path):
    nex_dic = {}
    with open(path) as nex_obj:
        for i in range(5):  # 跳过前五行
            nex_obj.readline()
        for line in nex_obj:
            # print(line)
            content = line.strip()
            if not content:
                break
            nex_dic[content.split('\t')[1]] = content.split('\t')[0]

    return nex_dic


# 解析分群文件 {样本名:群体信息}
def read_pop(path):
    pop_dic = {}
    with open(path) as pop_obj:
        for line in pop_obj:
            pop_dic[line.strip().split('\t')[0]] = line.strip().split('\t')[1]

# 解析表型数据 {样本名:表型值}
def reads_pheno(path):
    pheno_dic = {}
    with open(path) as pheno_obj:
        pheno_obj.readline()
        pheno_obj.tell()
        for line in pheno_obj:
            pheno_dic[line.strip().split('\t')[0]] = line.strip().split('\t')[1]
    return pheno_dic



# 工具函数，用于数据框匹配信息
def match(value,dic:dict):
    return dic.get(value,'NA')




# 主流程
def main(geno_file,nex_file,pheno_file,pop_file,out_file):
    out_header = ['Sample', 'Hap', 'Geno', 'Pheno', 'Population']
    # 读取SNP基因型文件
    snp_data = pd.read_csv(geno_file, sep='\t', header=None)
    snp_data.columns = ['Sample', 'Geno']

    # 解析单倍型文件、表型文件以及分群文件
    nex_dic = read_nex(nex_file)
    pheno_dic = reads_pheno(pheno_file)
    pop_dic = reads_pheno(pop_file)

    # 开始匹配
    snp_data['Hap'] = snp_data['Geno'].apply(match, args=(nex_dic,))
    snp_data['Pheno'] = snp_data['Sample'].apply(match, args=(pheno_dic,))
    snp_data['Population'] = snp_data['Sample'].apply(match, args=(pop_dic,))

    out_data = snp_data.loc[:, out_header]
    out_data.to_csv(out_file, index=False, header=True, sep='\t')


# 设置参数
def get_args():
    parser = argparse.ArgumentParser(description="转换vcf位点信息")
    parser.add_argument("-s", "--sample_geno", type=str,  help="样本名与基因型对应关系文件", required=True)
    parser.add_argument("-n", "--nex", type=str, help="单倍型名称与基因型对应关系文件", required=True)
    parser.add_argument("-phe", "--pheno", type=str, help="样本与表型对应关系文件", required=True)
    parser.add_argument("-pop", "--population", type=str, help="样本与分群对应关系文件", required=True)
    parser.add_argument("-o", "--out", type=str, help="输出文件名", required=True)
    return parser.parse_args()




if __name__ == '__main__':

    # 参数传递
    args = get_args()
    geno_file = args.sample_geno
    nex_file = args.nex
    pheno_file = args.pheno
    pop_file = args.population
    out_file = args.out
    main(geno_file,nex_file,pheno_file,pop_file,out_file)

    '''
    out_header = ['Sample','Hap','Geno','Pheno','Population']
    # 读取SNP基因型文件
    snp_data = pd.read_csv(geno_file,sep='\t',header=None)
    snp_data.columns = ['Sample','Geno']

    # 解析单倍型文件、表型文件以及分群文件
    nex_dic = read_nex(nex_file)
    pheno_dic = reads_pheno(pheno_file)
    pop_dic = reads_pheno(pop_file)

    # 开始匹配
    snp_data['Hap'] = snp_data['Geno'].apply(match,args=(nex_dic,))
    snp_data['Pheno'] = snp_data['Sample'].apply(match,args=(pheno_dic,))
    snp_data['Population'] = snp_data['Sample'].apply(match, args=(pop_dic,))

    out_data = snp_data.loc[:,out_header]
    out_data.to_csv(f'{out_file}.xls',index=False,header=True,sep='\t')
    '''





