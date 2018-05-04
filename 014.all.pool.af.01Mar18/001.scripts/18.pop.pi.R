### calculate pi within each population
### uses output from 12.calc.divergence.R
### making new script for this part for clarity
### DLF 24APR18

setwd("/lustre/scratch/projects/aquilegia/014.all.pool.af.01Mar18/002.input/")

### load divergence functions (also contains pi functions)
source("/lustre/scratch/projects/aquilegia/014.all.pool.af.01Mar18/001.scripts/00.divergence.functions.R")

load("12.calc.divergence.RData")
load("biallelic.allele.freq.bi.af.Rdata" )
pool.index <- read.table("pool.order.index.txt",stringsAsFactors=FALSE, header=TRUE)

##################################################
### get pi for each population  (by chr + genome)
###################################################

all.pop.pi <- matrix(NA, nrow(pool.index), ncol=8)
for (up.r in 1:nrow(pool.index)){
	print(up.r)
	n1 <- pool.index[up.r,6]*2
	up.div <- calc.pi(bi.af=bi.af, pop1=up.r, n1=n1, gsites=gsites, bi.snps=bi.snps)
	all.pop.pi[up.r,]<- up.div
	}

colnames(all.pop.pi)<- c(unique(bi.snps[,1]),"genome")
rownames(all.pop.pi) <- pool.index$pool.number
write.table(all.pop.pi, file="../003.output/all.pop.pi.txt")

### so, situation is pretty similar for all populations
### genome pi is a bit higher than in genomes paper (to be expected) - somewhere around 0.55-0.6
### chromosome four is eleveted in all
### chr3 is slightly lower in almost all (except GTT2)

###########################
### plot this by chromosome
############################## 

plot.var <- all.pop.pi[order(pool.index$species),] *100
pool.index.sort <- pool.index[order(pool.index$species),]
n.sp <- 10


pdf("../003.output/pop.pi.by.chr.pdf",width=7,height=5)
#quartz(width=7, height=5)
par(mai=c(1,1.5,0.5,0.75))
xlims=c(0,max(plot.var))
ylims=c(0.75,10.25)
chr.col <- c("blue","green","orange","red3","purple","aquamarine","violetred1")
plot(1:n.sp,xlim=xlims,ylim=ylims,type="l",xlab="percentage pairwise differences",yaxt="n",ylab="")
#rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "gray95")
abline(h=1:n.sp,col="gray80")
for (up in 1:n.sp){
	text(cbind(plot.var[up,1:7],up),labels=1:7,cex=1.5,col=chr.col,font=2)
}
axis(side=2, at=1:n.sp,labels=pool.index.sort$short.name,font=3,las=2)

### add species boxes
par(xpd=NA)
rect(0.91,0.55,0.95,5.45,col="red")
rect(0.91,5.55,0.95,10.45,col="lightgoldenrod")
text(0.93,3,labels="A. formosa",srt=90,adj=c(0.5,0.5),cex=1.1, font=3)
text(0.93,8,labels="A. pubescens",srt=90,adj=c(0.5,0.5),cex=1.1, font=3)
mtext("population",2,4.5, cex=1.2, font=2)
par(xpd=FALSE)

dev.off()
