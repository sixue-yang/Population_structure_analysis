# _author_ = 'sixueyang'
# _date_ = 2023/5/11 17:28
import sys
import os
import subprocess as sp
import pandas as pd
import argparse
from multiprocessing import Pool
import trun_vcf_v2
import match
import re
import threading
import time
import datetime
import plot_log

'''
input:
    gene_id
    gff
    总VCF
return:
    1、block与gene单倍型的箱线图
    2、block连锁强度 + 单倍型结构图
'''
Rscript = '/work1/Software/R3/bin/Rscript'
perl = '/work1/Users/wangtianyi/03.Software/01.Perl/perl/bin/perl'
plink = '/work1/Users/yangsixue/tools/conda/envs/bio/bin/plink'
nex_perl = '/work1/Users/jingxin/Pipeline/gwas_draw_gene_all_pip/getNEXV3.pl'
box_plot_R = '/work1/Users/yangsixue/pipline/GWAS_block_anlysis/script/boxplot.R'
LD_perl = '/work1/Users/yangsixue/pipline/GWAS_block_anlysis/script/gene_structure2.pl'


my_log = plot_log.Log()

def checkDir(path):
    if not os.path.exists(path):
        os.makedirs(path)

class Hap:
    def __init__(self,gff,gene_id,vcf_path,chrlist,pop_list,phe_path,outdir,pool:5,chunksize:10,step:5000,top:10):
        my_log.info('开始单倍型绘制')
        self.gff_path = gff
        self.chrlist = chrlist
        self.gene_id = gene_id
        self.vcf_path = vcf_path
        self.pop_list = pop_list
        self.phe_path = phe_path
        self.pool = pool
        self.chunksize = chunksize
        self.setp = step
        self.box_top = top
        self.filter_vcf_path = outdir + os.sep + '00.prepare' + os.sep + 'prepare.vcf.gz'
        self.gene_vcf_path = outdir + os.sep + '00.prepare' + os.sep + 'gene.vcf.gz'
        self.gene_vcf_pos_path = outdir + os.sep + '00.prepare' + os.sep + gene_id + 'pos'

        self.block_data = outdir + os.sep + '00.prepare' + os.sep + 'block' + os.sep + gene_id
        self.block_path = outdir + os.sep + '00.prepare' + os.sep + 'block' + os.sep + gene_id + '.blocks.det'
        self.block_vcf_path = outdir + os.sep + '00.prepare' + os.sep + 'block.vcf.gz'
        self.gene_pos = self.gene_vcf_pos_path + '.ldepth.mean'

        self.gene_result_dir = outdir + os.sep + '03.Hap' + os.sep + '01.gene'
        self.block_result_dir = outdir + os.sep + '03.Hap' + os.sep + '02.block'
        self.gene_nex = self.gene_result_dir + os.sep + 'gene.nex'
        self.block_nex = self.block_result_dir + os.sep + 'block.nex'
        self.gene_nex_list = self.gene_result_dir + os.sep + 'gene.snp.list'
        self.block_nex_list = self.block_result_dir + os.sep + 'block.snp.list'
        self.gene_boxplot_data = self.gene_result_dir + os.sep + 'gene.boxplot.xls'
        self.block_boxplot_data = self.block_result_dir + os.sep + 'block.boxplot.xls'
        self.gene_boxplot_pdf = self.gene_result_dir + os.sep + 'gene.boxplot.pdf'
        self.block_boxplot_pdf = self.block_result_dir + os.sep + 'block.boxplot.pdf'
        self.sample_list = outdir + os.sep + '00.prepare' + os.sep + 'sample.list'
        self.gene_ld_plot = self.gene_result_dir + os.sep + 'gene.LD.svg'
        self.block_ld_plot = self.block_result_dir + os.sep + 'block.LD.svg'
        checkDir(self.gene_result_dir)
        checkDir(self.block_result_dir)
        checkDir(os.path.dirname(self.gene_vcf_path))
        checkDir(os.path.dirname(self.block_data))
        self.mk_sample_list()
        self.chrs,self.star,self.end = self.get_gene_range()
        #print(self.chrs,self.star,self.end)
        self.chr_dic = self.trun2chr()




    # 提取特定染色体
    def vcf_filter(self):
        cmd = f'vcftools --gzvcf {self.vcf_path} --chr {self.chrs} --recode --recode-INFO-all -c | gzip -c > {self.filter_vcf_path}'
        my_log.info('提取vcf文件中特定染色体')
        sp.call(cmd,shell=True)

    def mk_sample_list(self):
        with open(self.pop_list) as pop_obj,open(self.sample_list, 'w') as sample_obj:
            for line in pop_obj:
                sample_obj.write(line.strip().split('\t')[0] + '\n')

    # 转换染色体转换表 -- > 字典
    def trun2chr(self):
        chr_dic = {}
        with open(self.chrlist) as chr_obj:
            for line in chr_obj:
                tag = line.strip().split('\t')
                chr_dic[tag[0]] = tag[1]
        return chr_dic

    # 获取基因所在位置
    def get_gene_range(self):
        chrs = 0
        star,end = 0,0
        with open(self.gff_path) as gff:
            for line in gff:
                if line.startswith("#"):
                    continue
                tag = line.strip().split('\t')
                types = tag[2]

                id_dic = dict(x.split('=') for x in tag[8].split(';'))
                # print(chrs,id_dic)
                if 'ID' in id_dic.keys():
                    if types == 'mRNA' and id_dic['ID'] == self.gene_id:
                        chrs = tag[0]
                        star, end = int(tag[3]),int(tag[4])
        if star:
            return chrs,star,end
        else:
            my_log.error('gene_id 不存在gff文件中，请检查gff文件格式')
            sys.exit('gene_id 不存在gff文件中，请检查gff文件格式')



    # 获取gene对应的vcf 并输出对应vcf的SNP位置信息
    def mk_gene_vcf(self):
        cmd = f'vcftools --gzvcf {self.filter_vcf_path} --chr {self.chrs} --from-bp {self.star} --to-bp {self.end} --recode --recode-INFO-all -c | gzip -c > {self.gene_vcf_path}'
        my_log.debug(f'获取gene区域对应的vcf:\n{cmd}')
        #print(cmd)
        sp.call(cmd,shell=True)
        pos_cmd = f'vcftools --gzvcf {self.gene_vcf_path} --site-mean-depth --out {self.gene_vcf_pos_path}'
        my_log.info(f'输出vcf的SNP位置信息:\n{pos_cmd}')
        sp.call(pos_cmd,shell=True)

    # 获取block信息
    def mk_block(self):
        cmd = f'sh /work1/Users/yangsixue/pipline/GWAS_block_anlysis/script/block.sh {self.filter_vcf_path} {self.block_data} {self.chrlist}'
        my_log.debug(f'计算block:\n{cmd}')
        sp.call(cmd,shell=True)

    # 根据基因范围选择block
    def get_merge_gene_block_pos(self):
        length = 0
        block_star,block_end = 0,0
        with open(self.block_path) as block_obj:

            gene_star = self.star - self.setp
            gene_end = self.end + self.setp
            for line in block_obj:

                if line.strip().startswith('CHR') or not line.strip():
                    continue
                else:

                    tag = re.split(r'\s+', line.strip())
                    # print(tag)
                    chrs = tag[5].split('|')[0].split(':')[0]
                    star,end = int(tag[1]),int(tag[2])

                    if chrs == self.chrs:
                        # print(star,end)
                        if star <= gene_star and end >= gene_end:
                            return star,end
                        elif star >= gene_star and end >= gene_end:
                            if gene_end - star > length:
                                length = gene_end - star
                                block_star,block_end = star,end
                        elif star <= gene_star and end <= gene_end:
                            if end - gene_star > length:
                                length = end - gene_star
                                block_star, block_end = star, end
                        elif star >= gene_star and end <= gene_end:
                            if end - star > length:
                                length = end - star
                                block_star, block_end = star, end
        return block_star, block_end

    # 提取block vcf
    def get_block_vcf(self,block_star,block_end):
        cmd = f'vcftools --gzvcf {self.filter_vcf_path} --chr {self.chrs} --from-bp {block_star} --to-bp {block_end} --recode --recode-INFO-all -c | gzip -c > {self.block_vcf_path}'
        my_log.debug(f'提取block区域对应vcf:\n{cmd}')
        sp.call(cmd,shell=True)

    # 单倍型分析
    def nex_analysis(self,vcf_path,out_prefix):
        cmd = f'perl {nex_perl} {vcf_path} {self.pop_list} 0 {out_prefix}'
        my_log.debug(f'计算单倍型:\n{cmd}')
        sp.call(cmd,shell=True)


    # gene_main
    def gene_main(self):
        # print(datetime.datetime.now(), '开始分析基因')
        my_log.info(f'开始分析基因单倍型')
        # 提取基因vcf
        my_log.info(f'提取基因vcf')
        self.mk_gene_vcf()
        # gene 单倍型分析
        my_log.info(f'gene 单倍型分析')
        self.nex_analysis(self.gene_vcf_path,self.gene_result_dir + os.sep + 'gene')
        # 转化vcf，生成单倍型
        my_log.info(f'转化基因vcf，生成单倍型')
        out_ped, out_map = trun_vcf_v2.mk_ped(self.gene_vcf_path, self.gene_result_dir + os.sep + 'gene', plink)
        trun_vcf_v2.read_ped(out_ped,self.chunksize,self.pool,self.gene_result_dir + os.sep + 'gene')

        # 匹配样本、单倍型、分群、表型
        my_log.info(f'匹配基因样本、单倍型、分群、表型')
        match.main(self.gene_nex_list,self.gene_nex,self.phe_path,self.pop_list,self.gene_boxplot_data)
        # 画箱线图

        cmd = f'{Rscript} {box_plot_R} --file {self.gene_boxplot_data} --out {self.gene_boxplot_pdf} --top {self.box_top}'
        # print(cmd)
        my_log.debug(f'绘制基因单倍型--表型箱线图:\n{cmd}')
        sp.call(cmd,shell=True)

        # 画连锁图
        ld_cmd = f'{perl} {LD_perl} {self.gff_path} {self.gene_id} {self.gene_vcf_path} {self.gene_result_dir + os.sep + "gene.rmfa"} {self.sample_list} > {self.gene_result_dir + os.sep + "gene.LD.svg"}'
        my_log.debug(f'绘制基因连锁图:\n{ld_cmd}')
        sp.call(ld_cmd,shell=True)
        my_log.info(f'基因绘图结束')
        # print(datetime.datetime.now(), '基因绘图结束')

    def block_main(self):
        # print(datetime.datetime.now(),'开始分析block')
        my_log.info(f'开始分析block单倍型')
        # 分析block 并提取block vcf
        my_log.info(f'开始分析block 并提取block vcf')
        self.mk_block()
        # 获取block所在区域
        block_star, block_end = self.get_merge_gene_block_pos()
        # print(block_star, block_end)
        if block_star:
            my_log.info(f'block区域存在，进行后续步骤。')
            # 获取block vcf文件
            my_log.info(f'开始block vcf文件。')
            self.get_block_vcf(block_star, block_end)
            # block 单倍型分析
            my_log.info(f'block单倍型分析。')
            self.nex_analysis(self.block_vcf_path, self.block_result_dir + os.sep + 'block')
            # 转化vcf --> 样本\t单倍型
            out_ped, out_map = trun_vcf_v2.mk_ped(self.block_vcf_path, self.block_result_dir + os.sep + 'block', plink)
            trun_vcf_v2.read_ped(out_ped, self.chunksize, self.pool, self.block_result_dir + os.sep + 'block')
            # 匹配样本、单倍型、分群、表型
            match.main(self.block_nex_list, self.block_nex, self.phe_path, self.pop_list, self.block_boxplot_data)
            # 画箱线图
            cmd = f'{Rscript} {box_plot_R} --file {self.block_boxplot_data} --out {self.block_boxplot_pdf} --top {self.box_top}'
            my_log.debug(f'block--单倍型箱线图绘制:\n{cmd}')
            sp.call(cmd,shell=True)

            ld_cmd = f'{perl} {LD_perl} {self.gff_path} {self.gene_id} {self.block_vcf_path} {self.block_result_dir + os.sep + "block.rmfa"} {self.sample_list} > {self.block_result_dir + os.sep + "block.LD.svg"}'
            my_log.debug(f'block连锁图绘制:\n{ld_cmd}')
            sp.call(ld_cmd,shell=True)
            my_log.info(f'block单倍型绘图结束')
            # print(datetime.datetime.now(), 'block绘图结束')
        else:
            my_log.warning(f'基因不在block区域内，请合理设置基因区范围')
            # print('基因不在block区域内，请合理设置基因区范围')



    # 主流程
    def main(self):
        # vcf 染色体过滤
        self.vcf_filter()
        t1 = threading.Thread(target=self.gene_main)
        t2 = threading.Thread(target=self.block_main)
        t1.start()
        t2.start()



if __name__ == '__main__':
    gff = sys.argv[1]
    gene_id = sys.argv[2]
    vcf_path = sys.argv[3]
    chrlist = sys.argv[4]
    pop_list = sys.argv[5]
    phe_path = sys.argv[6]
    outdir = sys.argv[7]

    run = Hap(gff,gene_id,vcf_path,chrlist,pop_list,phe_path,outdir,5,10,10000,10)
    run.main()










