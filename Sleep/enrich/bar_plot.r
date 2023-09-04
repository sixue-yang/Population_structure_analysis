#!/usr/bin/env Rscript
#@author: huangsonglin
#Email:huangsonglin@novogene.com
#jiangkai add graph_updown.png at 2019.1.30
suppressMessages({
library(reshape2)
library(ggplot2)
library(argparser)})

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
            newname <- paste(vectors[1],vectors[2],vectors[3],vectors[4],vectors[5],'...',sep=' ')
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

stat_frame <- read.delim(stat,head=TRUE,sep='\t',quote='')
up_down_stat<- stat
if (type=='align_pct') {
    stat <- stat_frame[,c('unique','multi','unmap')]
    stats <- data.frame()
    for (i in 1:nrow(stat)) {
        for (j in 1:ncol(stat)) {
            each_stat <- as.character(stat[i,j])
            stats[i,j] <- as.numeric(unlist(strsplit(unlist(strsplit(each_stat,'(',fixed=TRUE))[2],'%',fixed=TRUE))[1])
        }
    }
    colnames(stats) <- c('unique','multi','unmap')
    stat <- data.frame(sample=stat_frame$sample,stats)
    stat <- melt(stat)
    colnames(stat)<-c("sample","type",'percent')
    p <- ggplot(stat,aes(x=sample,y=percent,fill=type)) + geom_bar(stat='identity')
}

if (type=='align_region') {
    stat <- stat_frame[,c('exon','intron','intergenic')]
    stats <- data.frame()
    for (i in 1:nrow(stat)) {
        for (j in 1:ncol(stat)) {
            each_stat <- as.character(stat[i,j])
            stats[i,j] <- as.numeric(unlist(strsplit(unlist(strsplit(each_stat,'(',fixed=TRUE))[2],'%',fixed=TRUE))[1])
        }
    }
    colnames(stats) <- c('exon','intron','intergenic')
    stat <- data.frame(sample=stat_frame$sample,stats)
    stat <- melt(stat)
    colnames(stat)<-c("sample","type",'percent')
    p <- ggplot(stat,aes(x=sample,y=percent,fill=type)) + geom_bar(stat='identity')
}

if (type=='quant'|type=='biotype') {
    stat <- melt(stat_frame)
    colnames(stat)<-c("sample","type",'percent')
    p <- ggplot(stat,aes(x=sample,y=percent,fill=type)) + geom_bar(stat='identity')
}

if (type=='snp') {
    stat <- melt(stat_frame)
    colnames(stat)<-c("sample","type",'count')
    p <- ggplot(stat,aes(x=sample,y=count,fill=type)) + geom_bar(stat='identity')
}

if (type=='as') {
    stat <- melt(stat_frame)
    colnames(stat)<-c("compare","type",'ASevent')
    p <- ggplot(stat,aes(x=compare,y=ASevent,fill=type)) + geom_bar(stat='identity')
    if (nrow(stat_frame)==1) {
        p <- ggplot(stat,aes(x=compare,y=ASevent,fill=type)) + geom_bar(stat='identity',width=0.2)
    }
}

if (type=='das') {
    stat <- melt(stat_frame)
    colnames(stat)<-c("compare","type",'significant')
    p <- ggplot(stat,aes(x=compare,y=significant,fill=type)) + geom_bar(stat='identity')
    if (nrow(stat_frame)==1) {
        p <- ggplot(stat,aes(x=compare,y=significant,fill=type)) + geom_bar(stat='identity',width=0.2)
    }
}

if (type=='transcript') {
    stat <- melt(stat_frame)
    colnames(stat)<-c("sample","type",'transcript')
    p <- ggplot(stat,aes(x=sample,y=transcript,fill=sample)) + geom_bar(stat='identity')
}

if (type=='GO') {
    if("BP" %in% as.character(stat_frame$Category)){
    if (nrow(subset(stat_frame,Category=='BP'))<10){bprows=nrow(subset(stat_frame,Category=='BP'))}
    if (nrow(subset(stat_frame,Category=='BP'))>=10){bprows=10}}
    if("CC" %in% as.character(stat_frame$Category)){
    if (nrow(subset(stat_frame,Category=='CC'))<10){ccrows=nrow(subset(stat_frame,Category=='CC'))}
    if (nrow(subset(stat_frame,Category=='CC'))>=10){ccrows=10}}
    if("MF" %in% as.character(stat_frame$Category)){
    if (nrow(subset(stat_frame,Category=='MF'))<10){mfrows=nrow(subset(stat_frame,Category=='MF'))}
    if (nrow(subset(stat_frame,Category=='MF'))>=10){mfrows=10}}

    if (cutoff=='padj') {
        if("BP" %in% as.character(stat_frame$Category)){
        bp <- subset(stat_frame,Category=='BP')[1:bprows,c('Description','padj','Category','Count')]}else{bp <- data.frame()}
        if("CC" %in% as.character(stat_frame$Category)){
        cc <- subset(stat_frame,Category=='CC')[1:ccrows,c('Description','padj','Category','Count')]}else{cc <- data.frame()}
        if("MF" %in% as.character(stat_frame$Category)){
        mf <- subset(stat_frame,Category=='MF')[1:mfrows,c('Description','padj','Category','Count')]}else{mf <- data.frame()}
    } else if (cutoff=='pvalue') {
        if("BP" %in% as.character(stat_frame$Category)){
        bp <- subset(stat_frame,Category=='BP')[1:bprows,c('Description','pvalue','Category','Count')]}else{bp <- data.frame()}
        if("CC" %in% as.character(stat_frame$Category)){
        cc <- subset(stat_frame,Category=='CC')[1:ccrows,c('Description','pvalue','Category','Count')]}else{cc <- data.frame()}
        if("MF" %in% as.character(stat_frame$Category)){
        mf <- subset(stat_frame,Category=='MF')[1:mfrows,c('Description','pvalue','Category','Count')]}else{mf <- data.frame()}
    }

    stat <- rbind(bp,cc,mf)
    if(nrow(stat)==0){
        q()
    } else {
        Description <- rename(as.character(stat$Description))
        stat$Description <- factor(Description,levels=Description)

        if (cutoff=='padj') {
            p <- ggplot(stat,aes(x=Description,y=-log10(padj),fill=Category))
            p_bp <- ggplot(bp,aes(x=Description,y=-log10(padj),fill=Category))
            p_cc <- ggplot(cc,aes(x=Description,y=-log10(padj),fill=Category))
            p_mf <- ggplot(mf,aes(x=Description,y=-log10(padj),fill=Category))
        }
        if (cutoff=='pvalue') {
            p <- ggplot(stat,aes(x=Description,y=-log10(pvalue),fill=Category))
            p_bp <- ggplot(bp,aes(x=Description,y=-log10(pvalue),fill=Category))
            p_cc <- ggplot(cc,aes(x=Description,y=-log10(pvalue),fill=Category))
            p_mf <- ggplot(mf,aes(x=Description,y=-log10(pvalue),fill=Category))
        }
        p <- p + geom_bar(stat='identity') +theme(plot.margin=unit(c(1,1,2,4),'lines'))+geom_text(aes(label=Count,vjust=- 0.5),color="black",size=1.8)+theme(axis.text.x=element_text(hjust=1,angle=65,size=6.4,face ='bold'))
##jiangkai add
##jiangkai add graph_updown.png at 2019.1.30
        if (grepl('ALL',up_down_stat)){
            data <- read.delim(up_down_stat,header=T,sep="\t",quote='')
            str_trim <- function(str,num){
            if (nchar(str) > num){
                strTrim <- paste(substr(str,1,num),"...",sep="")
            } else {
                strTrim <- str
            }
            strTrim
        } 
        sub_data <- function(data,num){
	    if (nrow(data) > num){
	        x <- num
	    } else {
                x <- nrow(data)
	    }
            data[1:x,]
        } 

    if("BP" %in% as.character(data$Category)){
    BP_data2<-sub_data(subset(data,Category=="BP"),10)
    lables <- as.character(sapply(as.character(BP_data2[,3]),function(x){str_trim(x,35)}))
    BP_up_gene <- cbind(as.character(BP_data2[,1]),Ontology=lables,as.numeric(BP_data2[,11]),type=rep('up',nrow(BP_data2)))
    colnames(BP_up_gene) <- c('Ontology','Name','Number.of.Genes','type')
    BP_down_gene <- cbind(as.character(BP_data2[,1]),Ontology=lables,as.numeric(BP_data2[,13]),type=rep('down',nrow(BP_data2)))
    colnames(BP_down_gene) <- c('Ontology','Name','Number.of.Genes','type')}else{BP_up_gene<-c()
    BP_down_gene<-c()}

    if("CC" %in% as.character(data$Category)){
    CC_data2<-sub_data(subset(data,Category=="CC"),10)
    lables <- as.character(sapply(as.character(CC_data2[,3]),function(x){str_trim(x,35)}))
    CC_up_gene <- cbind(as.character(CC_data2[,1]),Ontology=lables,as.numeric(CC_data2[,11]),type=rep('up',nrow(CC_data2)))
    colnames(CC_up_gene) <- c('Ontology','Name','Number.of.Genes','type')
    CC_down_gene <- cbind(as.character(CC_data2[,1]),Ontology=lables,as.numeric(CC_data2[,13]),type=rep('down',nrow(CC_data2)))
    colnames(CC_down_gene) <- c('Ontology','Name','Number.of.Genes','type')}else{CC_up_gene<-c()
    CC_down_gene<-c()}

    if("MF" %in% as.character(data$Category)){
    MF_data2<-sub_data(subset(data,Category=="MF"),10)
    lables <- as.character(sapply(as.character(MF_data2[,3]),function(x){str_trim(x,35)}))
    MF_up_gene <- cbind(as.character(MF_data2[,1]),Ontology=lables,as.numeric(MF_data2[,11]),type=rep('up',nrow(MF_data2)))
    colnames(MF_up_gene) <- c('Ontology','Name','Number.of.Genes','type')
    MF_down_gene <- cbind(as.character(MF_data2[,1]),Ontology=lables,as.numeric(MF_data2[,13]),type=rep('down',nrow(MF_data2)))
    colnames(MF_down_gene) <- c('Ontology','Name','Number.of.Genes','type')}else{MF_up_gene<-c()
    MF_down_gene<-c()}

    BP_CC_MF_bind_data3=data.frame(rbind(BP_up_gene,BP_down_gene,CC_up_gene,CC_down_gene,MF_up_gene,MF_down_gene),stringAsFactors=F)
    result <- data.frame(BP_CC_MF_bind_data3[,1:2],Number.of.Genes=as.numeric(as.character(BP_CC_MF_bind_data3[,3])),type=BP_CC_MF_bind_data3[,4])

    P<-ggplot(result, aes(x=Name,fill=type,y=Number.of.Genes)) +  geom_bar(stat='identity',width=0.6 ,position="dodge")+facet_grid(.~Ontology,scale = "free")+ xlab("") + ylab(" Number of Genes ")+theme_bw()+theme(strip.text = element_text(size = rel(1.5)),strip.background=element_blank(),legend.text=element_text(size = rel(1)),legend.title = element_blank(),panel.grid.major=element_blank(),panel.grid.minor=element_blank(),panel.border=element_rect(colour = "black",size=1),panel.background=element_blank(),panel.grid =element_blank(),axis.line=element_line(),plot.title = element_text(size=20,face="bold",vjust = 1),axis.title.y=element_text(size =15,face ="bold",vjust = 0.3),axis.text.y=element_text(size =15, vjust = 1),axis.text.x= element_text(size =15, angle=80,face ="bold",vjust = 1, hjust =1))+scale_fill_manual(values = c("green","red"))
    ggsave(paste(prefix,'graph_updown.png',sep='_'),width=40,height=20,type='cairo-png',plot=P)
    ggsave(paste(prefix,'graph_updown.pdf',sep='_'),width=40,height=20,plot=P)
    }
}}

if (type=='KEGG'|type=='DO'|type=='Reactome'|type=='DisGeNET') {
    if(all(as.character(stat_frame[1,])=="NA")){
        q()
    } else {
        if (nrow(stat_frame)<20) {rows=nrow(stat_frame)}
        if (nrow(stat_frame)>=20) {rows=20}
        if (cutoff=='padj') {stat <- stat_frame[1:rows,c('Description','padj','Count')]}
        if (cutoff=='pvalue') {stat <- stat_frame[1:rows,c('Description','pvalue','Count')]}

        Description <- rename(as.character(stat$Description))
        stat$Description <- factor(Description,levels=Description)

        if (cutoff=='padj') {p <- ggplot(stat,aes(x=Description,y=-log10(padj)))}
        if (cutoff=='pvalue') {p <- ggplot(stat,aes(x=Description,y=-log10(pvalue)))}
        p <- p + geom_bar(stat='identity',fill='#6495ED') +
        theme(plot.margin=unit(c(1,1,2,4),'lines')) +
        theme(legend.position = "none")+geom_text(aes(label=Count,vjust=- 0.5),color="black",size=3)+theme(axis.text.x=element_text(hjust=1,angle=55,size=7.5,face ='bold'))
}}

p <- p + theme(panel.background=element_rect(fill="transparent"),axis.line=element_line())
#p <- p + theme(axis.text.x=element_text(hjust=1,angle=45,size=6))

pdf(file=paste(prefix,'.pdf',sep=''))
p
dev.off()

svg(filename=paste(prefix,'.svg',sep=''))
p
dev.off()
ggsave(file=paste(prefix,'.png',sep=''),type='cairo-png',plot=p)

if (type=='GO') {
    if("BP" %in% as.character(stat_frame$Category)){
    p_bp <- p_bp + geom_bar(stat='identity') +theme(plot.margin=unit(c(1,1,2,4),'lines'))+geom_text(aes(label=Count,vjust=- 0.5),color="black",size=1.8)+theme(axis.text.x=element_text(hjust=1,angle=65,size=6.4,face ='bold'))
    p_bp <- p_bp + theme(panel.background=element_rect(fill="transparent"),axis.line=element_line())
    #pdf(file=paste(prefix,'_bp.pdf',sep=''))
    #print(p_bp)
    #dev.off()
    svg(filename=paste(prefix,'_bp.svg',sep=''))
    print(p_bp)
    dev.off()
    ggsave(file=paste(prefix,'_bp.pdf',sep=''),plot=p_bp)
    ggsave(file=paste(prefix,'_bp.png',sep=''),type='cairo-png',plot=p_bp)
    #print ("HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH")
    }
     if("CC" %in% as.character(stat_frame$Category)){
    p_cc <- p_cc + geom_bar(stat='identity') +theme(plot.margin=unit(c(1,1,2,4),'lines'))+geom_text(aes(label=Count,vjust=- 0.5),color="black",size=1.8)+theme(axis.text.x=element_text(hjust=1,angle=65,size=6.4,face ='bold'))
    p_cc <- p_cc + theme(panel.background=element_rect(fill="transparent"),axis.line=element_line())
    #pdf(file=paste(prefix,'_cc.pdf',sep=''))
    #print(p_cc)
    #dev.off()
    svg(filename=paste(prefix,'_cc.svg',sep=''))
    print(p_cc)
    dev.off()
    ggsave(file=paste(prefix,'_cc.pdf',sep=''),plot=p_cc)
    ggsave(file=paste(prefix,'_cc.png',sep=''),type='cairo-png',plot=p_cc)
    }
    if("MF" %in% as.character(stat_frame$Category)){
    p_mf <- p_mf + geom_bar(stat='identity') +theme(plot.margin=unit(c(1,1,2,4),'lines'))+geom_text(aes(label=Count,vjust=- 0.5),color="black",size=1.8)+theme(axis.text.x=element_text(hjust=1,angle=65,size=6.4,face ='bold'))
    p_mf <- p_mf + theme(panel.background=element_rect(fill="transparent"),axis.line=element_line())
    #pdf(file=paste(prefix,'_mf.pdf',sep=''))
    #p_mf
    #dev.off()
    svg(filename=paste(prefix,'_mf.svg',sep=''))
    print(p_mf)
    dev.off()
    ggsave(file=paste(prefix,'_mf.pdf',sep=''),plot=p_mf)
    ggsave(file=paste(prefix,'_mf.png',sep=''),type='cairo-png',plot=p_mf)
    }
}
}
