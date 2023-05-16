# _author_ = 'sixueyang'
# _date_ = 2023/4/27 9:33
import argparse
import manhattanHap
import block
import datetime


if __name__ == '__main__':

    # 设置命令行参数解析器
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(dest='mode')
    # 创建子参数解析器1
    parser_all = subparsers.add_parser('all', help='绘制GWAS组合图(曼哈顿、倒三角、gene&block单倍型箱线图、连锁图)')
    parser_all.add_argument('-vcf', '--vcf', required=True, help='vcf文件')
    parser_all.add_argument('-m', '--mantattan', required=True, help='曼哈顿数据路径')
    parser_all.add_argument('-c', '--chrs', required=True, help='染色体编号')
    parser_all.add_argument('-s', '--snp_pos', required=True, help='受选择的 SNP 位置')
    parser_all.add_argument('-t', '--threshold', required=True, help='曼哈顿阈值线')
    parser_all.add_argument('-l', '--ranges', help='曼哈顿图 SNP 左右距离，默认10k',type=int, default=10000)
    parser_all.add_argument('-g', '--gff_path', required=True, help='基因结构注释文件 gff 路径')
    parser_all.add_argument('-gene', '--gene_id', required=True, help='关注的基因ID')
    parser_all.add_argument('-cl', '--chrlist', required=True, help='染色体转换表')
    parser_all.add_argument('-pop', '--pop_list', required=True, help='分群文件')
    parser_all.add_argument('-phe', '--phe_path', required=True, help='表型文件')
    parser_all.add_argument('-pool', '--pool', default=5,type=int, help='进程数量')
    parser_all.add_argument('-chunksize', '--chunksize', default=10, type=int,help='单倍型:一次解析的样本数量，默认一次10个样本进行解析。')
    parser_all.add_argument('-step', '--step', default=5000,type=int, help='单倍型:扩大基因区域，搜索基因所在block位置，针对挑选的基因不连锁的情况，默认扩大长度是上下5k。')
    parser_all.add_argument('-top', '--top', default=5, type=int,help='单倍型:绘制箱图时，频数20以上的单倍型全都存在杂合或缺失时，挑选频数最大数量的单倍型进行绘图，默认挑选top5')
    parser_all.add_argument('-o', '--outdir', required=True, help='输出目录')


    # 创建子参数解析器2
    parser_func1 = subparsers.add_parser('ManTriangles', help='绘制曼哈顿与倒三角')
    parser_func1.add_argument('-vcf', '--vcf', required=True, help='vcf文件')
    parser_func1.add_argument('-m', '--mantattan', required=True, help='曼哈顿数据路径')
    parser_func1.add_argument('-c', '--chrs', required=True, help='染色体编号')
    parser_func1.add_argument('-s', '--snp_pos', required=True, help='受选择的 SNP 位置')
    parser_func1.add_argument('-t', '--threshold', required=True, help='曼哈顿阈值线')
    parser_func1.add_argument('-l', '--ranges', help='曼哈顿图 SNP 左右距离，默认10k',default=10000)
    parser_func1.add_argument('-o', '--outdir', required=True, help='输出目录')

    # 创建子参数解析器3
    parser_func2 = subparsers.add_parser('plot_Hap', help='绘制基因&block所在单倍型箱图')
    parser_func2.add_argument('-vcf', '--vcf', required=True, help='vcf文件')
    parser_func2.add_argument('-g', '--gff_path', required=True, help='基因结构注释文件 gff 路径')
    parser_func2.add_argument('-gene', '--gene_id', required=True, help='关注的基因ID')
    parser_func2.add_argument('-cl', '--chrlist', required=True, help='染色体转换表')
    parser_func2.add_argument('-pop', '--pop_list', required=True, help='分群文件')
    parser_func2.add_argument('-phe', '--phe_path', required=True, help='表型文件')
    parser_func2.add_argument('-pool', '--pool', default=5, type=int, help='进程数量')
    parser_func2.add_argument('-chunksize', '--chunksize', default=10, type=int, help='单倍型:一次解析的样本数量，默认一次10个样本进行解析。')
    parser_func2.add_argument('-step', '--step', default=5000, type=int,
                            help='单倍型:扩大基因区域，搜索基因所在block位置，针对挑选的基因不连锁的情况，默认扩大长度是上下5k。')
    parser_func2.add_argument('-top', '--top', default=5, type=int,
                            help='单倍型:绘制箱图时，频数20以上的单倍型全都存在杂合或缺失时，挑选频数最大数量的单倍型进行绘图，默认挑选top5')
    parser_func2.add_argument('-o', '--outdir', required=True, help='输出目录')


    # 创建子参数解析器4
    # parser_func3 = subparsers.add_subparsers(dest='best')
    parser_func3 = subparsers.add_parser('one_Hap', help='单独绘制基因或block区域单倍型箱图')
    parser_func3.add_argument('-vcf', '--vcf', required=True, help='vcf文件')
    parser_func3.add_argument('-g', '--gff_path', required=True, help='基因结构注释文件 gff 路径')
    parser_func3.add_argument('-gene', '--gene_id', required=True, help='关注的基因ID')
    parser_func3.add_argument('-cl', '--chrlist', required=True, help='染色体转换表')
    parser_func3.add_argument('-pop', '--pop_list', required=True, help='分群文件')
    parser_func3.add_argument('-phe', '--phe_path', required=True, help='表型文件')
    parser_func3.add_argument('-pool', '--pool', default=5, type=int, help='进程数量')
    parser_func3.add_argument('-chunksize', '--chunksize', default=10, type=int, help='单倍型:一次解析的样本数量，默认一次10个样本进行解析。')
    parser_func3.add_argument('-step', '--step', default=5000, type=int,
                              help='单倍型:扩大基因区域，搜索基因所在block位置，针对挑选的基因不连锁的情况，默认扩大长度是上下5k。')
    parser_func3.add_argument('-top', '--top', default=5, type=int,
                              help='单倍型:绘制箱图时，频数20以上的单倍型全都存在杂合或缺失时，挑选频数最大数量的单倍型进行绘图，默认挑选top5')
    parser_func3.add_argument('-o', '--outdir', required=True, help='输出目录')
    parser_func3.add_argument('-type', '--type', required=True, help='选择绘制类型，默认只画gene [gene,block]',choices=['gene','block'],default='gene')


    args = parser.parse_args()

    # 根据`mode`参数选择要执行的操作
    if args.mode == 'all':
        vcf_path = args.vcf
        mantattan = args.mantattan
        chrs = args.chrs
        snp_pos = args.snp_pos
        threshold = args.threshold
        ranges = args.ranges
        outdir = args.outdir
        gff = args.gff_path
        gene_id = args.gene_id

        chrlist = args.chrlist
        pop_list = args.pop_list
        phe_path = args.phe_path
        pool = args.pool
        chunksize = args.chunksize
        step = args.step
        top = args.top
        print(datetime.datetime.now(), f"GWAS组合图开始绘制")
        man_run = manhattanHap.Mantattan(mantattan, chrs, snp_pos, ranges, threshold, outdir)
        man_run.plot_mantattan()
        print(datetime.datetime.now(), f"曼哈顿绘制完成")
        triangle_run = manhattanHap.Triangles(vcf_path, chrs, snp_pos, ranges, outdir)
        triangle_run.run()
        print(datetime.datetime.now(), f"倒三角绘制完成")
        hap_run = block.Hap(gff, gene_id, vcf_path, chrlist, pop_list, phe_path, outdir, pool, chunksize, step, top)
        hap_run.main()

        print(datetime.datetime.now(),f"所有图绘制完成")
    elif args.mode == 'ManTriangles':
        vcf_path = args.vcf
        mantattan = args.mantattan
        chrs = args.chrs
        snp_pos = args.snp_pos
        threshold = args.threshold
        ranges = args.ranges
        outdir = args.outdir
        print(datetime.datetime.now(), f"曼哈顿、倒三角开始绘制")
        run1 = manhattanHap.Mantattan(mantattan, chrs, snp_pos, ranges, threshold, outdir)
        run1.plot_mantattan()
        print(datetime.datetime.now(), f"曼哈顿绘制完成")
        run2 = manhattanHap.Triangles(vcf_path, chrs, snp_pos, ranges, outdir)
        run2.run()
        print(datetime.datetime.now(), f"倒三角绘制完成")

    elif args.mode == 'plot_Hap':

        vcf_path = args.vcf
        outdir = args.outdir
        gff = args.gff_path
        gene_id = args.gene_id

        chrlist = args.chrlist
        pop_list = args.pop_list
        phe_path = args.phe_path
        pool = args.pool
        chunksize = args.chunksize
        step = args.step
        top = args.top
        print(datetime.datetime.now(), f"基因与block单倍型相关图形开始绘制")
        hap_run = block.Hap(gff, gene_id, vcf_path, chrlist, pop_list, phe_path, outdir, pool, chunksize, step, top)
        hap_run.main()
        print(datetime.datetime.now(), f"绘制完成")
    elif args.mode == 'one_Hap':

        vcf_path = args.vcf
        outdir = args.outdir
        gff = args.gff_path
        gene_id = args.gene_id

        chrlist = args.chrlist
        pop_list = args.pop_list
        phe_path = args.phe_path
        pool = args.pool
        chunksize = args.chunksize
        step = args.step
        top = args.top
        type_choices = args.type
        hap_run = block.Hap(gff, gene_id, vcf_path, chrlist, pop_list, phe_path, outdir, pool, chunksize, step, top)
        hap_run.vcf_filter()
        if type_choices == 'gene':
            hap_run.gene_main()
        else:
            hap_run.block_main()
        print(datetime.datetime.now(),f" {type_choices} 绘制完成")