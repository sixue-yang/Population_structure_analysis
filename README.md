# Population_structure_analysis
用于记录群体结构分析过程中使用的小工具以及相关分析流程

## 目录说明
* tools 存放一些工具


## tools
merge_vcf.py 用于合并两个vcf文件，对两个vcf文件取SNP位置以及等位基因信息的交集

tips : bcftools同样能快速达到目的

```shell
bcftools isec A.vcf.gz B.vcf.gz -p dir
```

[bcftools参开文档]: https://www.jianshu.com/p/99895c7338b2



### GWAS_composite_graph_drawing/main.py

GWAS组合图绘制：曼哈顿图、倒三角、基因/block 对应的单倍型与表型箱线图、LD连锁图。

mian.py 设有不同的模式

**all**：绘制所有图

```python
options:
  -h, --help            show this help message and exit
  -vcf VCF, --vcf VCF   vcf文件
  -m MANTATTAN, --mantattan MANTATTAN
                        曼哈顿数据路径
  -c CHRS, --chrs CHRS  染色体编号
  -s SNP_POS, --snp_pos SNP_POS
                        受选择的 SNP 位置
  -t THRESHOLD, --threshold THRESHOLD
                        曼哈顿阈值线
  -l RANGES, --ranges RANGES
                        曼哈顿图 SNP 左右距离，默认10k
  -g GFF_PATH, --gff_path GFF_PATH
                        基因结构注释文件 gff 路径
  -gene GENE_ID, --gene_id GENE_ID
                        关注的基因ID
  -cl CHRLIST, --chrlist CHRLIST
                        染色体转换表
  -pop POP_LIST, --pop_list POP_LIST
                        分群文件
  -phe PHE_PATH, --phe_path PHE_PATH
                        表型文件
  -pool POOL, --pool POOL
                        进程数量
  -chunksize CHUNKSIZE, --chunksize CHUNKSIZE
                        单倍型:一次解析的样本数量，默认一次10个样本进行解析。
  -step STEP, --step STEP
                        单倍型:扩大基因区域，搜索基因所在block位置，针对挑选的基因不连锁的                                                                       情况，默认扩大长度是上下5k。
  -top TOP, --top TOP   单倍型:绘制箱图时，频数20以上的单倍型全都存在杂合或缺失时，挑选频数                                                                       最大数量的单倍型进行绘图，默认挑选top5
  -o OUTDIR, --outdir OUTDIR
                        输出目录
```

**ManTriangles**：仅绘制绘制曼哈顿与倒三角。

```python
usage: main.py ManTriangles [-h] -vcf VCF -m MANTATTAN -c CHRS -s SNP_POS -t THRESHOLD [-l RANGES] -o OUTDIR

options:
  -h, --help            show this help message and exit
  -vcf VCF, --vcf VCF   vcf文件
  -m MANTATTAN, --mantattan MANTATTAN
                        曼哈顿数据路径
  -c CHRS, --chrs CHRS  染色体编号
  -s SNP_POS, --snp_pos SNP_POS
                        受选择的 SNP 位置
  -t THRESHOLD, --threshold THRESHOLD
                        曼哈顿阈值线
  -l RANGES, --ranges RANGES
                        曼哈顿图 SNP 左右距离，默认10k
  -o OUTDIR, --outdir OUTDIR
                        输出目录

```

**plot_Hap**：仅绘制基因&block所在单倍型箱图以及各自对应的连锁

```python
usage: main.py plot_Hap [-h] -vcf VCF -g GFF_PATH -gene GENE_ID -cl CHRLIST -pop POP_LIST -phe PHE_PATH [-pool POOL] [-chunksize CHUNKSIZE] [-step STEP]
                        [-top TOP] -o OUTDIR

options:
  -h, --help            show this help message and exit
  -vcf VCF, --vcf VCF   vcf文件
  -g GFF_PATH, --gff_path GFF_PATH
                        基因结构注释文件 gff 路径
  -gene GENE_ID, --gene_id GENE_ID
                        关注的基因ID
  -cl CHRLIST, --chrlist CHRLIST
                        染色体转换表
  -pop POP_LIST, --pop_list POP_LIST
                        分群文件
  -phe PHE_PATH, --phe_path PHE_PATH
                        表型文件
  -pool POOL, --pool POOL
                        进程数量
  -chunksize CHUNKSIZE, --chunksize CHUNKSIZE
                        单倍型:一次解析的样本数量，默认一次10个样本进行解析。
  -step STEP, --step STEP
                        单倍型:扩大基因区域，搜索基因所在block位置，针对挑选的基因不连锁的情况，默认扩大长度是上下5k。
  -top TOP, --top TOP   单倍型:绘制箱图时，频数20以上的单倍型全都存在杂合或缺失时，挑选频数最大数量的单倍型进行绘图，默认挑选top5
  -o OUTDIR, --outdir OUTDIR
                        输出目录

```

**one_Hap**：单独绘制基因或block区域单倍型箱图以及各自对应的连锁

```python
usage: main.py one_Hap [-h] -vcf VCF -g GFF_PATH -gene GENE_ID -cl CHRLIST -pop POP_LIST -phe PHE_PATH [-pool POOL] [-chunksize CHUNKSIZE] [-step STEP]
                       [-top TOP] -o OUTDIR -type {gene,block}

options:
  -h, --help            show this help message and exit
  -vcf VCF, --vcf VCF   vcf文件
  -g GFF_PATH, --gff_path GFF_PATH
                        基因结构注释文件 gff 路径
  -gene GENE_ID, --gene_id GENE_ID
                        关注的基因ID
  -cl CHRLIST, --chrlist CHRLIST
                        染色体转换表
  -pop POP_LIST, --pop_list POP_LIST
                        分群文件
  -phe PHE_PATH, --phe_path PHE_PATH
                        表型文件
  -pool POOL, --pool POOL
                        进程数量
  -chunksize CHUNKSIZE, --chunksize CHUNKSIZE
                        单倍型:一次解析的样本数量，默认一次10个样本进行解析。
  -step STEP, --step STEP
                        单倍型:扩大基因区域，搜索基因所在block位置，针对挑选的基因不连锁的情况，默认扩大长度是上下5k。
  -top TOP, --top TOP   单倍型:绘制箱图时，频数20以上的单倍型全都存在杂合或缺失时，挑选频数最大数量的单倍型进行绘图，默认挑选top5
  -o OUTDIR, --outdir OUTDIR
                        输出目录
  -type {gene,block}, --type {gene,block}
                        选择绘制类型，默认只画gene [gene,block]
```

#### 快速测试

测试路径：/work1/Users/yangsixue/pipline/GWAS_block_anlysis/text/run_text/run.sh

```shell
python /work1/Users/yangsixue/pipline/GWAS_block_anlysis/script/main.py all -vcf /work1/Project/P20220125_ganlan_wangtianyi/07.gwas/SNP.filter.recode.vcf.gz -m /work1/Project/P20220125_ganlan_wangtianyi/07.gwas/snp/04.gwas_result/emmax/02.manhatan/IBE.mahattan.plot.data -c Chr8 -s 19097436 -t 5 -l 10000  -g /work1/Users/yangsixue/pipline/GWAS_block_anlysis/text/pipline_text/ysgl.HC.gff -gene C08p015830.1_BolYGL -cl /work1/Users/yangsixue/pipline/GWAS_block_anlysis/text/pipline_text/chr.list -pop /work1/Users/yangsixue/pipline/GWAS_block_anlysis/text/pipline_text/pop.list -phe text.phe -pool 5 -chunksize 10 -step 10000 -top 5 -o /work1/Users/yangsixue/pipline/GWAS_block_anlysis/text/run_text
```



## Sleep

群体项目中选择分析流程。

在原基础上进行环境及输出，XPCLR计算，KEGG、GO富集的优化。

======= 202401 优化内容 =======

* 染色体vcf进一步拆分，加速XPCLR计算
* KEGG、GO库，若直接输入参数，则不进行额外建库。若不存在，则单独建库
主脚本参数已修改。

============================



### 快速 运行

```shell
python sleep_main_v1.py -p [分群文件] -v [基础标记] -c [染色体编号] -b [基因的bed文件] -g [基因功能文件] -r [参考基因组] -S [KEGG物种缩写，默认为osa] -t [选择top的阈值,默认选取前百分之5,，0.05] -s [过滤SNP的阈值,默认过滤SNP数量小于2的窗口] -win [窗口长度] -step [步长] -job [最大任务数,默认100] -wt [qsub投递后等待结果时间,默认300秒] -cn [每条染色体进一步拆分的数量] -kl [KEGG数据库] -gl [go数据库] -sl [拆分数据tag标签输出路径]
```



#### 参数文件说明

##### 分群文件

```
# 第一列为样本，第二列为分群信息，列间使用tab键分隔
L266    P3
L267    P3
L297    P3
L189    P3
L905    P4
L1033   P2
L854    P2
L815    P2
L400    P2
L176    P2
L751    P2
L1050   P2
L1500   P7
L463    P7
L1495   P7
L51     P2
```

##### 基础标记

基本VCF文件，需要gz格式压缩

##### 染色体编号列表

```
# 染色体号
chr01
chr02
chr03
chr04
chr05
chr06
chr07
chr08
chr09
chr10
chr11
chr12
```

##### 基因的bed文件格式

```
# 染色体、起始位置、终止位置、基因，四列使用tab键分隔
chr01   2983    10815   Os01t0100100-01
chr01   11218   12435   Os01t0100200-01
chr01   11372   12284   Os01t0100300-00
chr01   12721   15685   Os01t0100400-01
chr01   12808   13978   Os01t0100466-00
chr01   16399   20144   Os01t0100500-01
chr01   22841   26892   Os01t0100600-01
chr01   25861   26424   Os01t0100650-00
chr01   27143   28644   Os01t0100700-01
chr01   29818   34453   Os01t0100800-01
chr01   35623   41136   Os01t0100900-01
chr01   58658   61090   Os01t0101150-00
chr01   60091   61086   Os01t0101175-00
chr01   62060   63576   Os01t0101200-01
chr01   62112   65537   Os01t0101200-02
```

##### 基因功能文件

```
# Gene_id Gene_id Function 列与列之间tab键分隔，表头如下
Gene_id Gene_id Function
Os01t0100100-01 Os01t0100100-01 sp|Q8BPQ7|SGSM1_MOUSE Small G protein signaling modulator 1 OS=Mus musculus GN=Sgsm1 PE=1 SV=2//7.06116e-22
Os01t0100200-01 Os01t0100200-01 sp|Q9LUC6|C7A14_ARATH Cytochrome P450 72A14 OS=Arabidopsis thaliana GN=CYP72A14 PE=2 SV=1//4.25192e-10
Os01t0100300-00 Os01t0100300-00 sp|Q9LUC6|C7A14_ARATH Cytochrome P450 72A14 OS=Arabidopsis thaliana GN=CYP72A14 PE=2 SV=1//4.12212e-10
Os01t0100400-01 Os01t0100400-01 sp|Q9SU40|SKU5_ARATH Monocopper oxidase-like protein SKU5 OS=Arabidopsis thaliana GN=SKU5 PE=1 SV=1//0
Os01t0100466-00 Os01t0100466-00 sp|Q9SU40|SKU5_ARATH Monocopper oxidase-like protein SKU5 OS=Arabidopsis thaliana GN=SKU5 PE=1 SV=1//6.98866e-16
Os01t0100500-01 Os01t0100500-01 -//-
Os01t0100600-01 Os01t0100600-01 -//-
Os01t0100650-00 Os01t0100650-00 -//-
```



#### 结果说明

```shell
01.Prepare/            # 各个分析模块分析前准备文件目录
02.Analysis/           # 各个分析模块分析结果目录
03.Top/                # Top 筛选后各个窗口的结果及基因
04.Enrich/             # 根据筛选后的基因进行KEGG及GO富集的结果
05.Result/             # 释放的结果目录
	[每一个群体组合]
		Arrange        # 存放几种方法的结果整理数据
			*total.value.xls  # 所有方法原始分析结果数据合并
			*top.value.xls    # 所有方法筛选后窗口数据合并结果
			*gene.stat.xls    # 基因维度，各种方法是否受选择及对应功能
			*win_sleep.xls    # 各个窗口，各种方法是否受选择统计
		Manhattan      # 各种选择方法曼哈顿图
		GO_KEGG        # 各种选择方法受选择基因富集结果
cluster_logs/  # qsub投递运行的日志输出

```





#### 参考测试路径

/work1/Users/yangsixue/pipline/group/sleep/snakemake3/run_test.sh







### 流程说明

本流程基于snakemake串写，具有断点续跑的功能，若对某一步骤结果不满意或需要调整，可以将那一部分目录删除或修改后，重新运行python主脚本。
