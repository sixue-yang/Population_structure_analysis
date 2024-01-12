# _author_ = 'sixueyang'
# _data_ = 2024/1/9 14:28
import gzip
import pandas as pd
import numpy as np
import os
import argparse

def split_vcf_gz(file_path, outdir,num_files=100):
    # Step 1: Count the header and content lines
    header_lines = []
    content_lines = 0
    with gzip.open(file_path, 'rt') as f:  # 'rt' mode for reading text from a gzip file
        for line in f:
            if line.startswith('#'):
                header_lines.append(line)
            else:
                content_lines += 1

    # Step 2: Calculate the number of lines per chunk
    lines_per_file = np.ceil(content_lines / num_files).astype(int)

    # Step 3: Read and split the content with pandas in chunks
    chunk_size = lines_per_file
    reader = pd.read_csv(file_path, compression='gzip', header=None, sep='\t', comment='#', chunksize=chunk_size)

    # Step 4: Write the header and chunks to new files
    out_prefix = outdir + '-split_'
    # log_list = []
    for i, chunk in enumerate(reader):

        # Create new VCF file for the chunk
        output_file = out_prefix + str(i + 1) + '.vcf.gz'
        with gzip.open(output_file, 'wt') as f:
            # Write the header
            f.writelines(header_lines)
            # Write the content
            chunk.to_csv(f, sep='\t', index=False, header=False,lineterminator='\n')

        # log_list.append(f'split_{i + 1}')

    # with open(log_path,'w') as log_obj:
    #     for i in log_list:
    #         print(i,file=log_obj)

def get_args():
    parser = argparse.ArgumentParser(description="选择分析 -- xpclr vcf文件拆分")
    parser.add_argument("-v", "--vcf", help="转换后vcf文件", required=True)
    parser.add_argument("-n", "--split_number", help="vcf文件拆分数量",type=int, default=100)
    # parser.add_argument("-l", "--log_path", help="日志路径", required=True)
    parser.add_argument("-o", "--out", help='输出文件所在目录', required=True)
    return parser.parse_args()


if __name__ == '__main__':
    args = get_args()

    file = args.vcf
    number_split = args.split_number
    # log_path = args.log_path
    outdir = args.out
    split_vcf_gz(file,outdir,number_split)



