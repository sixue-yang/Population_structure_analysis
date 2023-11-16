#!/usr/bin/env Rscript
#@author: huangsonglin
#Email:huangsonglin@novogene.com

suppressMessages({
library(argparser)
library(ggplot2)})

argv <- arg_parser('')
argv <- add_argument(argv,"--stat", help="the stat file")
argv <- add_argument(argv,"--type", help="the plot type")
argv <- add_argument(argv,"--cutoff", help="the cutoff value")
argv <- add_argument(argv,"--prefix", help="the prefix of outfile")
argv <- parse_args(argv)

stat <- argv$stat
type <- argv$type
cutoff <- argv$cutoff
prefix <- argv$prefix

rename <- function (x) {
    namelist <- c()
    for (i in 1:length(x)) {
        if (nchar(x[i]) >= 50) {
            vectors <- unlist(strsplit(x[i],' '))
            newname <- paste(vectors[1],vectors[2],vectors[3],vectors[4],'...',sep=' ')
            while(newname %in% namelist) {
                newname <- paste(newname,'.',sep='')
            }
            namelist[i] <- newname
        } else {
            namelist[i] <- x[i]
        }
    }
    return(namelist)
}

if(file.info(stat)$size == 0){}else{

stat <- read.delim(stat,head=TRUE,sep='\t',quote='',fill = TRUE )

if (type=='GO') {
    if("BP" %in% as.character(stat$Category)){
    if (nrow(subset(stat,Category=='BP'))<10){bprows=nrow(subset(stat,Category=='BP'))}
    if (nrow(subset(stat,Category=='BP'))>=10){bprows=10}}
    if("CC" %in% as.character(stat$Category)){
    if (nrow(subset(stat,Category=='CC'))<10){ccrows=nrow(subset(stat,Category=='CC'))}
    if (nrow(subset(stat,Category=='CC'))>=10){ccrows=10}}
    if("MF" %in% as.character(stat$Category)){
    if (nrow(subset(stat,Category=='MF'))<10){mfrows=nrow(subset(stat,Category=='MF'))}
    if (nrow(subset(stat,Category=='MF'))>=10){mfrows=10}}

    if (cutoff=='padj') {
    if("BP" %in% as.character(stat$Category)){
    bp <- subset(stat,Category=='BP')[1:bprows,c('GeneRatio','Description','padj','Count')]}else{bp <- data.frame()}
    if("CC" %in% as.character(stat$Category)){
    cc <- subset(stat,Category=='CC')[1:ccrows,c('GeneRatio','Description','padj','Count')]}else{cc <- data.frame()}
    if("MF" %in% as.character(stat$Category)){
    mf <- subset(stat,Category=='MF')[1:mfrows,c('GeneRatio','Description','padj','Count')]}else{mf <- data.frame()}
    } else if (cutoff=='pvalue') {
    if("BP" %in% as.character(stat$Category)){
    bp <- subset(stat,Category=='BP')[1:bprows,c('GeneRatio','Description','pvalue','Count')]}else{bp <- data.frame()}
    if("CC" %in% as.character(stat$Category)){
    cc <- subset(stat,Category=='CC')[1:ccrows,c('GeneRatio','Description','pvalue','Count')]}else{cc <- data.frame()}
    if("MF" %in% as.character(stat$Category)){
    mf <- subset(stat,Category=='MF')[1:mfrows,c('GeneRatio','Description','pvalue','Count')]}else{mf <- data.frame()}
    }

    stat <- rbind(bp,cc,mf)
    if(nrow(stat)!=0){
        ratio <- matrix(as.numeric(unlist(strsplit(as.character(stat$GeneRatio),"/"))),ncol=2,byrow=TRUE)
        stat$GeneRatio <- ratio[,1]/ratio[,2]

        Description <- rename(as.character(stat$Description))
        stat$Description <- factor(Description,levels=Description)

        if (cutoff=='padj') {
           stat <- stat[order(stat$padj,decreasing=F),]
            p <- ggplot(stat,aes(x=GeneRatio,y=Description,colour=padj,size=Count))
        } else if (cutoff=='pvalue') {
            stat <- stat[order(stat$pvalue,decreasing=F),]
           p <- ggplot(stat,aes(x=GeneRatio,y=Description,colour=pvalue,size=Count))
        }

    p <- p + geom_point() +
        scale_colour_gradientn(colours=rainbow(4),guide="colourbar") + 
        expand_limits(color=seq(0,1,by=0.25)) +
        theme(axis.text=element_text(color="black",size=10))
   } else {
     q()
   }
}

if (type=='KEGG'|type=='DO'|type=='Reactome'|type=='DisGeNET') {
    if(all(as.character(stat[1,])=="NA")){
       q() 
    } else {
        if (nrow(stat)<20) {rows=nrow(stat)}
        if (nrow(stat)>=20) {rows=20}
        if (cutoff=='padj') {stat <- stat[1:rows,c('GeneRatio','Description','padj','Count')]}
        if (cutoff=='pvalue') {stat <- stat[1:rows,c('GeneRatio','Description','pvalue','Count')]}

        ratio <- matrix(as.numeric(unlist(strsplit(as.character(stat$GeneRatio),"/"))),ncol=2,byrow=TRUE)
        stat$GeneRatio <- ratio[,1]/ratio[,2]

        Description <- rename(as.character(stat$Description))
        stat$Description <- factor(Description,levels=Description)

        if (cutoff=='padj') {
            stat <- stat[order(stat$padj,decreasing=F),]
            p <- ggplot(stat,aes(x=GeneRatio,y=Description,colour=padj,size=Count))
        } else if (cutoff=='pvalue') {
            stat <- stat[order(stat$pvalue,decreasing=F),]
            p <- ggplot(stat,aes(x=GeneRatio,y=Description,colour=pvalue,size=Count))
        }

        p <- p + geom_point() +
        scale_colour_gradientn(colours=rainbow(4),guide="colourbar") +
        expand_limits(color=seq(0,1,by=0.25)) +
        theme(axis.text=element_text(color="black",size=10))
    }
}

p <- p + theme_bw() + theme(panel.border=element_rect(colour="black"))

ggsave(p,filename=paste(prefix,'.pdf',sep=''))


svg(filename=paste(prefix,'.svg',sep=''))
p
dev.off()
ggsave(filename=paste(prefix,'.png',sep=''),type='cairo-png',plot=p)
}
