# _author_ = 'sixueyang'
# _date_ = 2023/8/28 22:47
from Bio import SeqIO
from Bio.Seq import Seq
import argparse
from collections import defaultdict



def parse_gene_bed(path):
    gene_dic = defaultdict(lambda: defaultdict(list))
    with open(path) as bed_obj:
        for line in bed_obj:
            chrs,start,end,gene_id = line.strip().split('\t')
            # 过滤异常基因 cds 碱基数不是3的倍数
            if (int(end) - int(start) + 1) % 3 == 0:
                gene_dic[chrs][gene_id] = [int(start)-1,int(end)]
            else:
                continue
    return gene_dic


# 设置参数
def get_args():
    parser = argparse.ArgumentParser(description="选择分析 -- 富集 转换CDS与pep序列.")
    parser.add_argument("-f", "--fasta", help="参考基因组", required=True)
    parser.add_argument("-b", "--bed", help=r'包含基因位置信息的bed文件，无表头 chrs\tStar\tEnd\tGeneid', required=True)
    parser.add_argument("-oc", "--out_cds", help='输出cds文件路径', required=True)
    parser.add_argument("-op", "--out_pep", help='输出pep文件路径', required=True)
    return parser.parse_args()


# 处理蛋白序列 替换末尾的*号 将中间的终止密码子替换成X
def trun_pep(seq:str):
    new_seq = seq.rstrip('*')
    new_seq = new_seq.replace('*','X')
    return new_seq


if __name__ == '__main__':

    argv = get_args()
    fasta_file = argv.fasta
    gene_bed_file = argv.bed
    out_cds = argv.out_cds
    out_pep = argv.out_pep


    gene_dic = parse_gene_bed(gene_bed_file)
    # 创建一个结果列表，用于存储提取出的蛋白序列
    with open(fasta_file) as fasta_file_obj,open(out_cds,'w') as out_cds_obj,open(out_pep,'w') as out_pep_obj:

        # 遍历 FASTA 文件中的每个记录
        for record in SeqIO.parse(fasta_file, "fasta"):
            gene_sequence = record.seq
            for gene in gene_dic[str(record.id)]:

                gene_start,gene_end = gene_dic[str(record.id)][gene]  # 设置基因起始位置
                # 提取基因序列
                gene_seq = gene_sequence[gene_start:gene_end]
                # 将基因序列转化为蛋白序列
                protein_seq = Seq(str(gene_seq)).translate()
                protein_seq = trun_pep(protein_seq)
                print(f'>{gene}\n{gene_seq}',file=out_cds_obj,flush=True)
                print(f'>{gene}\n{protein_seq}',file=out_pep_obj,flush=True)





