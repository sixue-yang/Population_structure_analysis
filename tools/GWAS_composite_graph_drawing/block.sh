in=$1
prefix=$2
chr=$3
perl /work1/Users/jingxin/Pipeline/vcf2ped/vcf2pedMap.pl $in $chr $prefix
/work1/Software/Plink/bin/plink  --file $prefix --blocks no-pheno-req --blocks-max-kb 1000 --blocks-min-maf 0.05 --blocks-strong-lowci 0.70 --blocks-strong-highci 0.98 --blocks-recomb-highci 0.90 --blocks-inform-frac 0.95 --out $prefix
