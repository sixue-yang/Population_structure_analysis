args=commandArgs(T)
p1=read.table(args[1])
pos=p1$V2/1000000
val=p1$V3
start=as.numeric(args[2])
end=as.numeric(args[3])
pos_line=as.numeric(args[6])
xlim<-c(start/1000000,end/1000000)
ymax<-max(val) + 1
pdf(args[7])
#png(paste0(args[7], ".png"), width = 900, height = 600)
plot(x=pos,y=val,xlim=xlim,ylim=c(0,ymax),xlab=args[4],ylab='-log10(P)', pch=16,col = "blue",cex=0.5,cex.axis=1.0, cex.lab=1.0,cex.main=1.0)
abline(h=args[5],type="b",lty=2,col="grey",xpd=F,lwd=1.5)
abline(v=pos_line/1000000,type="b",lty=2,col="grey",xpd=F,lwd=1.5)
dev.off()

png_filename = sub("\\.pdf$", ".png", args[7])
png(png_filename, width = 1600, height = 1300)
plot(x=pos,y=val,xlim=xlim,ylim=c(0,ymax),xlab=args[4],ylab='-log10(P)', pch=16,col = "blue",cex=0.5,cex.axis=1.0, cex.lab=1.0,cex.main=1.0)
abline(h=args[5],type="b",lty=2,col="grey",xpd=F,lwd=1.5)
abline(v=pos_line/1000000,type="b",lty=2,col="grey",xpd=F,lwd=1.5)

dev.off()
