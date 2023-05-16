# _author_ = 'sixueyang'
# _date_ = 2023/4/27 9:51
import pandas as pd
import os
import subprocess as sp
import sys

# 环境固定参数
Rscript = '/work1/Software/R3/bin/Rscript'
vcftools = '/work1/Software/vcftools-1.15/bin/vcftools'
perl = '/work1/Users/wangtianyi/03.Software/01.Perl/perl/bin/perl'
haploview = '/work1/Software/haploview/bin/haploview'
plink = '/work1/Software/Plink/bin/plink'
# 脚本固定参数设置
trans_perl = '/work1/Users/jingxin/Pipeline/gwas_haplotype_mahateen_draw/trans2.pl'
# 曼哈顿图绘制
Mantattan_plot_R = '/work1/Users/jingxin/Pipeline/gwas_draw_gene_all_pip/draw.region.manhateen.R'
# block 倒三角
Block_triangle_plot_R = '/work1/Users/yangsixue/pipline/GWAS_block_anlysis/script/LDheatmap.R'


# 文件上级目录判断

def checkDir(path):
    if not os.path.exists(path):
        os.makedirs(path)


'''
功能:绘制SNP区域的曼哈顿图
input:
    mantattan: GWAS分析曼哈顿绘图结果路径
    chrs: 染色体位置
    snp_pos: SNP位置
    range: SNP前后范围 若想整个窗口的长度为8kb，则需要输入4000000
    threshold: 阈值
    outfile: 输出目录
return:
    目录:01.HapMantattan
    曼哈顿图
'''

class Mantattan:
    def __init__(self,mantattan,chrs,snp_pos,range,threshold,outdir):
        self.mantattan = mantattan
        self.chrs = chrs
        self.snp_pos = snp_pos
        self.star = int(snp_pos) - int(range)
        self.end = int(snp_pos) + int(range)
        self.threshold = threshold
        self.out_dir = outdir + os.sep + '01.Mantattan'
        self.outmantattan = self.out_dir + os.sep + 'mantattan.pdf'
        checkDir(self.out_dir)

        self.draw_path = self.out_dir + os.sep + "mantattan.draw.data"
        _df = self.data_filter()
        _df.to_csv(self.draw_path,sep='\t',index=False,header=None)

    # 曼哈顿数据依据染色体号及区域过滤
    def data_filter(self):
        data = pd.read_csv(self.mantattan,sep='\t',header=None)
        return data[(data[0]==self.chrs)&(data[1]>=self.star)&(data[1]<=self.end)]
    # 绘图
    def plot_mantattan(self):
        cmd = f'{Rscript} {Mantattan_plot_R} {self.draw_path} {self.star} {self.end} {self.chrs} {self.threshold} {self.snp_pos} {self.outmantattan}'
        sp.call(cmd,shell=True)


'''
功能:绘制SNP区域的倒三角
input:
    vcfpath: GWAS分析曼哈顿绘图结果路径
    chrs: 染色体位置
    snp_pos: SNP位置
    range: SNP前后范围 若想整个窗口的长度为8kb，则需要输入4000000
return:
    目录:01.HapMantattan
    倒三角
'''

def trun_num(value):
    if value == '-9':
        return '0'
    elif value == 'A':
        return '1'
    elif value == 'C':
        return '2'
    elif value == 'G':
        return '3'
    elif value == 'T':
        return '4'
    else:
        return value

class Triangles:
    def __init__(self,vcfpath,chrs,snp_pos,range,outdir):
        self.vcf_path = vcfpath
        self.chrs = chrs
        self.snp_pos = snp_pos
        self.star = int(snp_pos) - int(range)
        self.end = int(snp_pos) + int(range)
        self.out_dir = outdir + os.sep + '02.Triangles'
        checkDir(self.out_dir)
        self.out_file = self.out_dir + os.sep + 'triangles'


    def run(self):
        # vcftools 提取SNP位点并输出为ped与map
        vcf2ped_cmd = f'{vcftools} --gzvcf {self.vcf_path} --chr {self.chrs} --from-bp {self.star} --to-bp {self.end} --plink-tped --out {self.out_file}'
        ped_trun_cmd = f'{plink} --tfile {self.out_file} --recode --out {self.out_file}'

        sp.call(vcf2ped_cmd,shell=True)
        sp.call(ped_trun_cmd, shell=True)
        with open(f'{self.out_file}.map') as map_file,open(f'{self.out_file}.info','w') as info,open(f'{self.out_file}.ped') as ped_obj,open(f'{self.out_file}.trans.ped','w') as trans_obj:
            # 转换map至info
            for line in map_file:
                pos = line.strip().split('\t')[-1]
                info.write(f'1:{pos}\t{pos}\n')
            # 转换ped文件中ATCG碱基值
            for line in ped_obj:
                content = ' '.join(map(trun_num, line.strip().split(' ')))
                print(content, file=trans_obj)
        haploview_cmd = f'''{haploview} -memory 20000 -maxdistance 30000 -n -log {self.out_file}.log -pedfile {self.out_file}.trans.ped -hwcutoff 0 -out {self.out_file} -dprime   -minMAF  0.01 -png'''
        trans1_cmd = f'cut -f 1 {self.out_file}.LD |sort -u|sort -t ":" -k1,1 -k2n,2  >{self.out_file}.LD.mark'
        trans2_cmd = f'{perl} {trans_perl} {self.out_file}.LD {self.out_file}.LD.mark {self.out_file}.LD.input'
        triangles_plot = f'{Rscript} {Block_triangle_plot_R} --file {self.out_file}.LD.input --out {self.out_file}.block_LDheatmap.png'
        print(triangles_plot)
        cmd_list = [haploview_cmd,trans1_cmd,trans2_cmd,triangles_plot]
        for cmds in cmd_list:
            sp.call(cmds,shell=True)



if __name__ == '__main__':

    vcf_path = sys.argv[1]
    mantattan = sys.argv[2]
    chrs = sys.argv[3]
    snp_pos = sys.argv[4]
    threshold = sys.argv[5]
    ranges = sys.argv[6]
    outdir = sys.argv[7]
    run1 = Mantattan(mantattan,chrs,snp_pos,ranges,threshold,outdir)
    run1.plot_mantattan()


    run2 = Triangles(vcf_path,chrs,snp_pos,ranges,outdir)
    run2.run()






