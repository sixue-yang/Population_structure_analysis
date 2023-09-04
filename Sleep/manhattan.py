#############################
##zhangyazhou 2023.08.16
#############################
import sys

import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import re,os
import datetime
import argparse

datetime_object = datetime.datetime.now()
parser = argparse.ArgumentParser()

parser.add_argument('-i', "--input", metavar="", required=True, help="输入曼哈顿输入文件")
parser.add_argument('-ts', "--threshold", metavar="", default="", help="输入阈值线")
parser.add_argument('-label_s', "--label_s", metavar="", type=float, default=12,
                    help="输入图标题大小")
parser.add_argument('-sc', "--sign_color", metavar="", default="",
                    help="输入qq图和曼哈顿图超过阈值线颜色，\"#FF0000\"")
parser.add_argument('-mc', "--manhattan_color", metavar="", default="#0F6AB1,#D86814,#25823B",
                    help="输入曼哈顿颜色，\"#f07400,#0088f0\"")
parser.add_argument('-msw', "--manhattan_scatter_kw", metavar="", type=float, default=12,
                    help="输入曼哈顿图点大小,默认是12")
parser.add_argument('-mtw', "--manhattan_threshold_w", metavar="", type=float, default=1,
                    help="输入曼哈顿图阈值线宽度,默认是1")
parser.add_argument('-mtc', "--manhattan_threshold_c", metavar="", default="r",
                    help="输入曼哈顿图点大小，\"#FF0000\"")
parser.add_argument('-pw', "--plt_w", metavar="", type=float, default=12, help="输出图片宽度,单位是英寸")
parser.add_argument('-ph', "--plt_h", metavar="", type=float, default=4, help="输出图片高度,单位是英寸")
parser.add_argument('-o', "--output", metavar="", default="./", help="输出文件，默认输出./")
args = parser.parse_args()

os.makedirs(args.output, exist_ok=True)
def deal_color(color):
    if "\"" in color:
        color = color.eval()
    else:
        color = color
    return color


# 取染色体的数字列
def get_new_lis(item):
    new_list = []
    # 使用正则表达式提取数字部分
    for i in item:
        number = re.findall(r'\d+', i)
        if number:
            new_list.append(str(int("".join(number))))
        else:
            new_list.append(i)  # 如果没有数字，则将其原来的
    return new_list


#字体设置
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42
# 正常显示字体，和正常显示负号
plt.rcParams['font.sans-serif'] = ['FZHei-B01S']
plt.rcParams['axes.unicode_minus'] = False


######################################绘制曼哈顿图################################
def manhattan(data, chr_list, sign_color, threshold, manhattan_threshold_w, manhattan_threshold_c, manhattan_scatter_kw,
              label_s, name, CHRS):
    # 设置图片长宽
    fig = plt.figure(figsize=(args.plt_w, args.plt_h), dpi=300)
    #print("""#run 绘制曼哈顿图""")
    color_set = args.manhattan_color.strip().split(",")
    color = []
    for i in range(len(chr_list)):
        res = (i + 1) % len(color_set)
        color.append(color_set[res])
    # 开始做图:
    # 处理x轴上每个chr标注的位置，这里使用每条染色体上SNP索引的中间值
    chrom_df = data.groupby(CHRS)['i'].median()
    # 在上一步计算得到的位置添加x轴的刻度和标注
    plot = sns.scatterplot(data=data, x='i', y=name, hue=CHRS, palette=color, linewidth=0, s=manhattan_scatter_kw,
                           legend=None)
    if sign_color:
        new_data = data[data[name].gt(threshold)]
        new_color = []
        for i in range(len(new_data.drop_duplicates(subset=[CHRS])[CHRS].tolist())):
            new_color.append(deal_color(sign_color))
        sns.scatterplot(data=new_data, x='i', y=name, hue=CHRS, palette=new_color, linewidth=0,
                        s=manhattan_scatter_kw, legend=None)
    plot.set_xticks(chrom_df)
    if len(chr_list) > 12:
        plot.set_xticklabels(get_new_lis(chrom_df.index), fontsize=label_s, rotation=0)
    else:
        plot.set_xticklabels(chrom_df.index, fontsize=label_s, rotation=0)
    plt.yticks(fontsize=10)
    plot.spines['left'].set_linewidth(1)
    plot.spines['bottom'].set_linewidth(1)
    # 最后添加各种标注，作为阈值的判断线，调整图像Y轴范围
    plot.set_xlabel('Chromosome', fontsize=label_s, labelpad=6)
    plot.set_ylabel(name, fontsize=label_s, labelpad=6)
    plot.spines["right"].set_visible(False)
    plot.spines["top"].set_visible(False)
    # plt.title(args.prefix + " manhattan plot", fontsize=title_s)
    plot.axhline(y=threshold, linewidth=manhattan_threshold_w, linestyle="--",
                 color=deal_color(manhattan_threshold_c))
    plot.margins(0.02, 0)
    plot.set_ylim(0, data[name].max() * 1.2)
    plt.savefig(args.output + "/" + name + ".png", bbox_inches='tight', pad_inches=0.3)
    plt.savefig(args.output + "/" + name + ".pdf", bbox_inches='tight', pad_inches=0.3)


def deal_data():
    #读取数据
    datas = pd.read_csv(args.input, header=0, sep="\t")
    #读取分析的名字
    name_list = datas.columns.values.tolist()
    name_lis = name_list[4:]
    datas['sort'] = list(map(int, get_new_lis(datas[name_list[0]].tolist())))
    datas = datas.sort_values(by=['sort', name_list[1]])
    #替换inf数据
    datas.replace([np.inf, -np.inf], np.nan, inplace=True)
    #datas = datas.reset_index(drop=True)
    if len(name_lis) == 0:
        sys.exit("您输入的文件数据不正确")
    #计算位置
    datas['pos'] = datas[[name_list[1], name_list[2]]].mean(axis=1).round().astype('int')
    #循环画图
    for name in name_lis:
        new_data = datas[[name_list[0], "pos", name]].copy()
        #new_data.mask(new_data[name]< 0, 0)
        new_data.dropna(inplace=True)
        new_data = new_data.reset_index(drop=True)
        new_data['i'] = new_data.index
        chr_list = new_data.drop_duplicates(subset=[name_list[0]])[name_list[0]].tolist()
        if not args.threshold:
            if len(new_data[name].sort_values(ascending=False).tolist())>1:
                threshold = new_data[name].sort_values(ascending=False).tolist()[int(len(new_data) * 0.05) - 1]
            else:
                threshold = 0
        else:
            threshold = float(args.threshold)
        if len(new_data[name].sort_values(ascending=False).tolist())>1:
            manhattan(new_data, chr_list, args.sign_color, threshold, args.manhattan_threshold_w,
                  args.manhattan_threshold_c, args.manhattan_scatter_kw, args.label_s, name, name_list[0])


def man():
    deal_data()


if __name__ == '__main__':
    man()

