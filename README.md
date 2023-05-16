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



​	
