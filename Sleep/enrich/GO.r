#!/usr/bin/env Rscript
#@author: huangsonglin
#Email:huangsonglin@novogene.com
#jiangkai add cnet plot at 2018.12.2

suppressMessages({
library(clusterProfiler)
library(argparser)})

argv <- arg_parser('')
argv <- add_argument(argv,"--go", help="the go annotation file")
argv <- add_argument(argv,"--diffgene", help="the differential gene")
argv <- add_argument(argv,"--diffresult", help="the differential analysis result")
argv <- add_argument(argv,"--enrich", help="the enrich method")
argv <- add_argument(argv,"--prefix", help="the prefix of outfile")
argv <- parse_args(argv)

go <- argv$go
diffgene <- argv$diffgene
diffresult <- argv$diffresult
enrich <- argv$enrich
prefix <- argv$prefix
compare_list <-unlist(strsplit(prefix,split='/'))
compare<-compare_list[length(compare_list)]

go_df <- read.delim(go,header=TRUE,sep='\t',quote='')
term2gene_bp <- go_df[go_df[,'go_ontology']=='BP',c('go_id','gene_id')]
term2name_bp <- go_df[go_df[,'go_ontology']=='BP',c('go_id','go_term')]
term2gene_cc <- go_df[go_df[,'go_ontology']=='CC',c('go_id','gene_id')]
term2name_cc <- go_df[go_df[,'go_ontology']=='CC',c('go_id','go_term')]
term2gene_mf <- go_df[go_df[,'go_ontology']=='MF',c('go_id','gene_id')]
term2name_mf <- go_df[go_df[,'go_ontology']=='MF',c('go_id','go_term')]

if (enrich=='normal') {
    diffgene_df <- read.delim(diffgene,header=TRUE,sep='\t',quote='')
    diffgene_vt <- as.character(diffgene_df$Gene)
    genename_vt <- as.character(diffgene_df$Function)
    names(genename_vt) <- diffgene_vt
    if (!is.na(diffresult)) {
        diffresult_df <- read.delim(diffresult,header=TRUE,sep='\t',quote='')
        backgene_vt <- as.character(diffresult_df$Gene)
        #backgene_foldchage <- as.numeric(diffresult_df$log2FoldChange)
        #names(backgene_foldchage) <- backgene_vt
        BPenrich <- enricher(gene=diffgene_vt,universe=backgene_vt,TERM2GENE=term2gene_bp,TERM2NAME=term2name_bp,pAdjustMethod='BH',pvalueCutoff=1,qvalueCutoff=1)
        BPenrich@ontology="BP"
        pdf(paste(prefix,".GObp_DAG.pdf",sep=""))
        try(plotGOgraph(BPenrich,firstSigNodes = 5))
        dev.off()
        if (file.info(paste(prefix,".GObp_DAG.pdf",sep=""))$size==3611) {file.remove(paste(prefix,".GObp_DAG.pdf",sep=""))}
        png(paste(prefix,".GObp_DAG.png",sep=""),type="cairo-png",width=480*4,height=480*4,res=72*4)
        try(plotGOgraph(BPenrich,firstSigNodes = 5))
        dev.off()

        CCenrich <- enricher(gene=diffgene_vt,universe=backgene_vt,TERM2GENE=term2gene_cc,TERM2NAME=term2name_cc,pAdjustMethod='BH',pvalueCutoff=1,qvalueCutoff=1)
        CCenrich@ontology="CC"
        pdf(paste(prefix,".GOcc_DAG.pdf",sep=""))
        try(plotGOgraph(CCenrich,firstSigNodes = 5))
        dev.off()
        if (file.info(paste(prefix,".GOcc_DAG.pdf",sep=""))$size==3611) {file.remove(paste(prefix,".GOcc_DAG.pdf",sep=""))}
        png(paste(prefix,".GOcc_DAG.png",sep=""),type="cairo-png",width=480*4,height=480*4,res=72*4)
        try(plotGOgraph(CCenrich,firstSigNodes = 5))
        dev.off()

        MFenrich <- enricher(gene=diffgene_vt,universe=backgene_vt,TERM2GENE=term2gene_mf,TERM2NAME=term2name_mf,pAdjustMethod='BH',pvalueCutoff=1,qvalueCutoff=1)
        MFenrich@ontology="MF"
        pdf(paste(prefix,".GOmf_DAG.pdf",sep=""))
        try(plotGOgraph(MFenrich,firstSigNodes = 5))
        dev.off()
        if (file.info(paste(prefix,".GOmf_DAG.pdf",sep=""))$size==3611) {file.remove(paste(prefix,".GOmf_DAG.pdf",sep=""))}
        png(paste(prefix,".GOmf_DAG.png",sep=""),type="cairo-png",width=480*4,height=480*4,res=72*4)
        try(plotGOgraph(MFenrich,firstSigNodes = 5))
        dev.off()

    } else {
        BPenrich <- enricher(gene=diffgene_vt,TERM2GENE=term2gene_bp,TERM2NAME=term2name_bp,pAdjustMethod='BH',pvalueCutoff=1,qvalueCutoff=1)
        CCenrich <- enricher(gene=diffgene_vt,TERM2GENE=term2gene_cc,TERM2NAME=term2name_cc,pAdjustMethod='BH',pvalueCutoff=1,qvalueCutoff=1)
        MFenrich <- enricher(gene=diffgene_vt,TERM2GENE=term2gene_mf,TERM2NAME=term2name_mf,pAdjustMethod='BH',pvalueCutoff=1,qvalueCutoff=1)
    }

    BPenrich <- as.data.frame(BPenrich)
    if (length(BPenrich) >1){
    names(BPenrich)[6] <- 'padj'
    BPenrich <- subset(BPenrich,select=-qvalue)
    Category <- rep('BP',time=nrow(BPenrich))
    BPenrich <- cbind(Category,BPenrich)
    }
    CCenrich <- as.data.frame(CCenrich)
    if (length(CCenrich) >1){
    names(CCenrich)[6] <- 'padj'
    CCenrich <- subset(CCenrich,select=-qvalue)
    Category <- rep('CC',time=nrow(CCenrich))
    CCenrich <- cbind(Category,CCenrich)
    }
    MFenrich <- as.data.frame(MFenrich)
    if (length(MFenrich) >1){
    names(MFenrich)[6] <- 'padj'
    MFenrich <- subset(MFenrich,select=-qvalue)
    Category <- rep('MF',time=nrow(MFenrich))
    MFenrich <- cbind(Category,MFenrich)
    }
    GOenrich <- rbind(BPenrich,CCenrich,MFenrich)
    if (length(GOenrich) >1){
    names(GOenrich)[2] <- 'GOID'
    gene_list <- strsplit(GOenrich$geneID,split='/')
    geneName <- unlist(lapply(gene_list,FUN=function(x){paste(genename_vt[x],collapse='/')}))
    GOenrich <- data.frame(GOenrich[,1:8],geneName,GOenrich[,9,drop=F])
    Significant <- subset(GOenrich,padj<0.05)
    }else{
    Significant <- NULL
    }
} else if (enrich=='GSEA') {
    diffresult_df <- read.delim(diffresult,header=TRUE,sep='\t',quote='')
    genename_vt <- as.character(diffresult_df$Gene)
    names(genename_vt) <- as.character(diffresult_df$Function)
    foldchange <- diffresult_df[,'log2FoldChange']
    names(foldchange) <- as.character(diffresult_df$Gene)
    foldchange <- sort(foldchange,decreasing=TRUE)

    BPenrich <- GSEA(geneList=foldchange,TERM2GENE=term2gene_bp,TERM2NAME=term2name_bp,pvalueCutoff=1,pAdjustMethod="BH")
    CCenrich <- GSEA(geneList=foldchange,TERM2GENE=term2gene_cc,TERM2NAME=term2name_cc,pvalueCutoff=1,pAdjustMethod="BH")
    MFenrich <- GSEA(geneList=foldchange,TERM2GENE=term2gene_mf,TERM2NAME=term2name_mf,pvalueCutoff=1,pAdjustMethod="BH")

    BPenrich <- as.data.frame(BPenrich)
    names(BPenrich)[7] <- 'padj'
    BPenrich <- subset(BPenrich,select=-qvalues)
    Category <- rep('BP',time=nrow(BPenrich))
    BPenrich <- cbind(Category,BPenrich)

    CCenrich <- as.data.frame(CCenrich)
    names(CCenrich)[7] <- 'padj'
    CCenrich <- subset(CCenrich,select=-qvalues)
    Category <- rep('CC',time=nrow(CCenrich))
    CCenrich <- cbind(Category,CCenrich)

    MFenrich <- as.data.frame(MFenrich)
    names(MFenrich)[7] <- 'padj'
    MFenrich <- subset(MFenrich,select=-qvalues)
    Category <- rep('MF',time=nrow(MFenrich))
    MFenrich <- cbind(Category,MFenrich)

    GOenrich <- rbind(BPenrich,CCenrich,MFenrich)
    names(GOenrich)[2] <- 'GOID'
    gene_list <- strsplit(GOenrich$core_enrichment,split='/')
    geneName <- unlist(lapply(gene_list,FUN=function(x){paste(genename_vt[x],collapse='/')}))
    GOenrich <- data.frame(GOenrich,geneName)
    Significant <- subset(GOenrich,padj<0.05)
}

write.table(GOenrich,file=paste(prefix,'_GOenrich.xls',sep=''),sep='\t',quote=F,row.names=F)
write.table(Significant,file=paste(prefix,'_GOenrich_significant.xls',sep=''),sep='\t',quote=F,row.names=F)

#GOenrich.xls add up down gene
if (grepl('_deg_all.xls',diffgene)){
    split=unlist(strsplit(diffgene,split="_"))
    combine=paste(split[1:(length(split)-2)],collapse="_")
    if(length(unlist(strsplit(compare,split="vs")))==2){
    diffgene_up <- as.character(read.delim(paste(combine,'_deg_up.xls',sep=""),header=T,sep="\t")$Gene,quote='')
    diffgene_down <- as.character(read.delim(paste(combine,'_deg_down.xls',sep=""),header=T,sep="\t")$Gene,quote='')
    GOenrich=read.delim(paste(prefix,'_GOenrich.xls',sep=''),header=T,quote='')
    geneID <- strsplit(as.character(GOenrich$geneID),split="/")
    inter_geneup_id_num=unlist(lapply(geneID,function(x){length(intersect(x,diffgene_up))}))
    inter_geneup_id=lapply(geneID,function(x){intersect(x,diffgene_up)})
    inter_geneup_id_combine <- unlist(lapply(inter_geneup_id,function(x){paste(x,collapse='/')}))
    inter_genedown_id_num=unlist(lapply(geneID,function(x){length(intersect(x,diffgene_down))}))
    inter_genedown_id=lapply(geneID,function(x){intersect(x,diffgene_down)})
    inter_genedown_id_combine <- unlist(lapply(inter_genedown_id,function(x){paste(x,collapse='/')}))
    out=data.frame(GOenrich,Up=inter_geneup_id_num,Up_Gene_id=inter_geneup_id_combine,Down=inter_genedown_id_num,Down_Gene_id=inter_genedown_id_combine)
    write.table(out,file=paste(prefix,'_GOenrich.xls',sep=''),sep='\t',quote=F,row.names=F)
#GOenrich_significant.xls add up down gene

    GOenrich=read.delim(paste(prefix,'_GOenrich_significant.xls',sep=''),header=T,quote='')
    geneID <- strsplit(as.character(GOenrich$geneID),split="/")
    inter_geneup_id_num=unlist(lapply(geneID,function(x){length(intersect(x,diffgene_up))}))
    inter_geneup_id=lapply(geneID,function(x){intersect(x,diffgene_up)})
    inter_geneup_id_combine <- unlist(lapply(inter_geneup_id,function(x){paste(x,collapse='/')}))
    inter_genedown_id_num=unlist(lapply(geneID,function(x){length(intersect(x,diffgene_down))}))
    inter_genedown_id=lapply(geneID,function(x){intersect(x,diffgene_down)})
    inter_genedown_id_combine <- unlist(lapply(inter_genedown_id,function(x){paste(x,collapse='/')}))
    out=data.frame(GOenrich,Up=inter_geneup_id_num,Up_Gene_id=inter_geneup_id_combine,Down=inter_genedown_id_num,Down_Gene_id=inter_genedown_id_combine)
    write.table(out,file=paste(prefix,'_GOenrich_significant.xls',sep=''),sep='\t',quote=F,row.names=F)
}}
