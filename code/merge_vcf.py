# _author_ = 'sixueyang'
# _date_ = 2023/4/23 11:21
import vcf
import sys
import time
import os
import subprocess as sp
import argparse
from collections import defaultdict
import threading
import multiprocessing
from multiprocessing import Process
from queue import Queue

m = multiprocessing.Manager()
result_queue = m.Queue()

# vcf 数据前处理
class Prepare:

    def __init__(self,vcf_path):
        self.vcf_path = vcf_path
        print(vcf_path)
        self.dirs = os.path.dirname(self.vcf_path)
        self.filter_name = os.path.basename(self.vcf_path).replace('vcf.gz', 'filter.vcf.gz')
        self.out_path = os.path.join(self.dirs, self.filter_name)
        self.file_name = os.path.basename(vcf_path).split('.')[0] + '.txt'


    # vcf 基础过滤 (先过滤，节约后续数据处理时间) 返回输出文件路径
    def filter_vcf(self):
        sample_num = self.get_samples_num()
        if sample_num > 500:
            maf = 0.01
        else:
            maf = 0.05

        cmd = f'{vcftools} --gzvcf {self.vcf_path} --max-missing 0.8 --maf {maf} --remove-indels --max-alleles 2 --min-alleles 2 --minQ 20 --recode --recode-INFO-all -c | gzip -c > {self.out_path}'
        sp.call(cmd,shell=True)



    # 获取vcf样本数量
    def get_samples_num(self):
        vcf_obj = vcf.Reader(filename=self.vcf_path)
        return len(vcf_obj.samples)


    # 生成vcf文件前五列信息文件 返回输出位置信息路径
    def mk_pos(self):


        if self.filter_name.endswith('gz'):
            cmd = f'zcat {self.filter_name} | cut -f 1-5 > {self.file_name}'
        else:
            cmd = f'cut -f 1-5 {self.filter_name} > {self.file_name}'
        sp.call(cmd,shell=True)

    # 解析位置信息
    def parse_pos(self):
        pos_set = []
        with open(self.file_name) as file_obj:
            for line in file_obj:
                if line.startswith('#'):
                    continue
                CHROM, POS, ID, REF, ALT = line.strip().split('\t')
                pos_set.append(','.join(map(str, [CHROM, POS, REF, ALT])))
        # print(pos_set)
        result_queue.put(pos_set)

        # print(result_queue.qsize())



    def main(self):
        # vcf基础过滤
        self.filter_vcf()
        # 截取前五列
        self.mk_pos()
        # 解析位置信息
        self.parse_pos()





# 将vcf数据存入字典
def get_snp(vcf_path,pos_list:list,dic:dict):
    vcf_obj = vcf.Reader(filename=vcf_path)

    for record in vcf_obj:
        if f'{record.CHROM}:{record.POS}' in pos_list:
            dic[f'{record.CHROM}:{record.POS}'].append(record)
    return dic




# 合并
def merge_record(vcf1_path,vcf2_path,out_name,data_dic:dict):
    vcf1_obj = vcf.Reader(filename=vcf1_path)
    vcf2_obj = vcf.Reader(filename=vcf2_path)
    vcf1_obj.samples = vcf1_obj.samples + vcf2_obj.samples

    out_vcf = vcf.Writer(open(f'{out_name}.vcf.gz','w'),vcf1_obj)
    # out_vcf.write_header(in_vcf._header)
    for chr_pos in data_dic.keys():
        record1,record2 = data_dic[chr_pos]
        # 构建sample_indexs
        sample_indexs = {}
        i = 0
        for sample in record1.samples + record1.samples:
            sample_indexs[sample.sample] = i
            i += 1


        merged_record = vcf.model._Record(
            record1.CHROM,
            record1.POS,
            record1.ID,
            record1.REF,
            record1.ALT,
            record1.QUAL,
            record1.FILTER,
            record1.INFO,
            record1.FORMAT,
            sample_indexs,
            samples=record1.samples + record2.samples
        )
        out_vcf.write_record(merged_record)
    out_vcf.close()



class MyProcess(Process): #继承Process类
    def __init__(self,file):
        super(MyProcess,self).__init__()
        self.file = file


    def run(self):
        vcf_obj = Prepare(self.file)
        vcf_obj.main()


# 设置参数开关
def get_args():
    parser = argparse.ArgumentParser(description="合并两个vcf文件，取交集")
    parser.add_argument("-f1", "--vcf1", help="第一个vcf路径", required=True)
    parser.add_argument("-f2", "--vcf2", help="第二个vcf路径", required=True)
    parser.add_argument("-o", "--out_name", help="输出文件头名", required=True)
    parser.add_argument("-vcftools", "--vcftools", help="vcftools路径", default='/work1/Users/yangsixue/tools/conda/bin/vcftools')
    return parser.parse_args()



if __name__ == '__main__':

    args = get_args()
    # 输入的两个vcf文件
    file1 = args.vcf1
    file2 = args.vcf2
    out_name = args.out_name
    vcftools = args.vcftools

    file_list = [file1,file2]
    # 两个文件位置信息列表
    set_list = []
    # 两个进程进行数据前处理
    process_list = []

    for i,file in enumerate(file_list):
        p = MyProcess(file)  # 实例化进程对象
        p.start()
        # print(f'正在进行第{i}个vcf文件预处理...')
        process_list.append(p)

    for i in process_list:
        i.join()
        # set_list.append(result_queue.get())
    set_list = [set(result_queue.get()) for i in file_list]
    # print(set_list)
    print(f'完成预处理工作')
    unio_set = set_list[0] & set_list[1]
    unio_list = [f'{chr_pos.split(",")[0]}:{chr_pos.split(",")[1]}' for chr_pos in unio_set]
    # 输出交集、测试用
    text_obj = open('pos.list', 'w')
    for i in unio_list:
        text_obj.write(i + '\n')

    print('开始vcf数据提取')
    # 定义数据字典
    out_dic = defaultdict(list)

    # 先处理文件1、再处理文件2 顺序改变会影响后续合并写入时的顺序颠倒
    out_dic = get_snp(file1, unio_list, out_dic)
    out_dic = get_snp(file2, unio_list, out_dic)
    print('开始合并vcf数据')
    # 合并文件
    merge_record(file1, file2, out_name, out_dic)
    print('合并完毕')
    # print(len(vcf_obj.samples))
    # print('成功读取vcf句柄',time.time() - times)




