# _author_ = 'sixueyang'
# _date_ = 2023/8/14 16:44
import argparse
import gzip

# TODO:转换vcf格式，输入文件默认gz格式，若有其它格式需求，请联系研发部
# 工具函数
def trun(value:str):
    return value.split(':')[0]


# 设置参数
def get_args():
    parser = argparse.ArgumentParser(description="选择分析 -- xpclr软件输入vcf文件格式转换.")
    parser.add_argument("-v", "--vcf", help="原始vcf文件", required=True)
    parser.add_argument("-o", "--out", help='输出文件路径', required=True)
    return parser.parse_args()


# 主函数
def main():
    argv = get_args()
    vcf_path = argv.vcf
    out_file = argv.out

    with gzip.open(vcf_path,'rt') as inf, gzip.open(out_file,'wt') as outf:
        for line in inf:
            if line.startswith('#'):
                outf.write(line)
                outf.flush()
            else:
                tag = line.strip().split('\t')
                tag[8:] = list(map(trun,tag[8:]))
                # print(tag)
                outf.write('\t'.join(tag)+'\n')
                outf.flush()


if __name__ == '__main__':
    main()


