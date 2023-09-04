# _author_ = 'sixueyang'
# _date_ = 2023/9/1 14:03
import sys
import os
import argparse
import yaml
import subprocess as sp



# 创建qsub投递日志目录
def check_out_Dir(yaml_path):
    with open(yaml_path, 'r') as config_file:
        config_data = yaml.safe_load(config_file)
        for rule in config_data:
            if 'output' in config_data[rule]:
                strout = os.path.join(os.path.dirname(os.path.dirname(config_data[rule]['output'])), rule)
                os.makedirs(strout,exist_ok=True)


# 设置传参
def get_args():
    parser = argparse.ArgumentParser(description="选择分析 - 主控制流程")
    parser.add_argument("-p", "--pop_file", help="群体分群文件", required=True)
    parser.add_argument("-v", "--vcf_file", help="vcf文件 gz压缩格式", required=True)
    parser.add_argument("-c", "--chr_list", help="染色体列表",required=True)
    parser.add_argument("-b", "--bed_file", help="包含基因位置信息的bed文件",required=True)
    parser.add_argument("-g", "--gene_name_file", help="基因对应的功能注释文件", required=True)
    parser.add_argument("-r", "--ref_fasta", help="参考基因组", required=True)
    parser.add_argument("-S", "--spe", help="KEGG 建库使用的物种缩写", default='osa')
    parser.add_argument("-t", "--topCut", help="筛选top的百分比，默认筛选前百分之五", type=float,default=0.05)
    parser.add_argument("-s", "--SNPCut", help="过滤无效窗口的SNP阈值，默认过滤2以下", type=int,default=2)
    parser.add_argument("-win", "--win", help="窗口大小 bp", type=int,required=True)
    parser.add_argument("-step", "--step", help="步长 bp", type=int,required=True)
    parser.add_argument("-job", "--max_jobs", help="最大并行任务数，默认100", type=int,default=100)
    parser.add_argument("-wt", "--wait_time", help="任务执行后等待结果的时间，默认300秒", type=int, default=300)
    parser.add_argument("-sn", "--snakemake", help="snakemake软件路径", default='/work1/Users/yangsixue/tools/conda/envs/snakemake/bin/snakemake')
    return parser.parse_args()

if __name__ == '__main__':
    argv = get_args()
    pop_file = argv.pop_file
    vcf_file = argv.vcf_file
    chr_list = argv.chr_list
    bed_file = argv.bed_file
    gene_name_file = argv.gene_name_file
    ref_fasta = argv.ref_fasta
    spe = argv.spe
    topCut = argv.topCut
    SNPCut = argv.SNPCut
    win = argv.win
    step = argv.step
    job = argv.max_jobs
    wait_time = argv.wait_time
    snakemake = argv.snakemake

    cluster_contig_path = os.path.join(os.path.dirname(sys.argv[0]),'cluster.yaml')
    # 创建投递标准输出目录
    check_out_Dir(cluster_contig_path)
    # 写项目配置文件
    input_configfile = {
        'pop_file': pop_file,
        'vcf_file': vcf_file,
        'chr_list': chr_list,
        'bed_file': bed_file,
        'gene_name_file': gene_name_file,
        'fast_file': ref_fasta,
        'spe': spe,
        'topCut': topCut,
        'SNPCut': SNPCut,
        'win': win,
        'step': step,
        'script_dir':os.path.dirname(sys.argv[0])
    }
    with open('input_config.yaml', 'w') as config_file:
        yaml.dump(input_configfile, config_file)

    snakemake_script = os.path.join(os.path.dirname(sys.argv[0]),'Sleep.snakemake')
    sub_cmd = f'''{snakemake} -s {snakemake_script} --configfile input_config.yaml --cluster-config {cluster_contig_path} --cluster "qsub -V -cwd -l vf={{cluster.memory}},p={{cluster.threads}} -o {{cluster.output}} -e {{cluster.error}}" --jobs {job} --latency-wait {wait_time}
    '''
    print('命令',sub_cmd)
    sp.call(sub_cmd,shell=True)
