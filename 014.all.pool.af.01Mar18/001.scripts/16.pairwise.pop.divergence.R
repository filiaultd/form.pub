### calculate divergence between all population pairs
### uses output from 12.calc.divergence.R
### making new script for this part for clarity
### DLF 23APR18



setwd("/lustre/scratch/projects/aquilegia/014.all.pool.af.01Mar18/002.input/")

### load divergence functions
source("/lustre/scratch/projects/aquilegia/014.all.pool.af.01Mar18/001.scripts/00.divergence.functions.R")

load("12.calc.divergence.RData")  ### this is from old name of script 12.
load("biallelic.allele.freq.bi.af.Rdata" )
pool.index <- read.table("pool.order.index.txt",stringsAsFactors=FALSE, header=TRUE)

#### get comparisons to make
up.comps <- expand.grid(1:ncol(bi.af),1:ncol(bi.af))
up.comps <- up.comps[up.comps[,1]!=up.comps[,2],]
#### still some duplicates here!  check this!
up.co <- apply(up.comps,1, function(x){
	x <- x[order(x)]
})
up.co <- t(up.co)
up.co <- unique(up.co)

### get divergence for each comparison (by chr + genome)

all.pair.div <- matrix(NA, nrow(up.co), ncol=8)
for (up.r in 1:nrow(up.co)){
	up.dat <- up.co[up.r,]
	up.dat <- as.numeric(up.dat[order(up.dat)])
	print(up.dat)
	pop1 <- pool.index[up.dat[1],3]
	pop2 <- pool.index[up.dat[2],3]
	n1 <- pool.index[up.dat[1],6]*2
	n2 <- pool.index[up.dat[2],6]*2	
	up.div <- calc.div(bi.af=bi.af, pop1=pop1, pop2=pop2, n1=n1, n2=n2, gsites=gsites, bi.snps=bi.snps)
	all.pair.div[up.r,]<- up.div
	}

colnames(all.pair.div)<- c(unique(bi.snps[,1]),"genome")

write.table(all.pair.div, file="../003.output/all.pairwise.divergence.txt")
#all.pair.div <- read.table(file="../003.output/all.pairwise.divergence.txt",stringsAsFactors=FALSE)

################################	
# between vs within divergence?
##################################

colnames(up.co) <- c("Var1","Var2")
comp.div <- cbind(up.co, all.pair.div) 
comp.div <- merge(comp.div, pool.index[,3:5], by.x="Var1", by.y="order") 
colnames(comp.div)[11:12]<-paste("Var1",colnames(comp.div)[11:12], sep=".")
comp.div <- merge(comp.div, pool.index[,3:5], by.x="Var2", by.y="order") 
colnames(comp.div)[13:14]<-paste("Var2",colnames(comp.div)[13:14], sep=".")
comp.div$type <- with(comp.div, paste(Var1.species, Var2.species, sep="."))
for (up in 1:nrow(comp.div)){
	if (comp.div[up,12]==comp.div[up,14]){comp.div[up,15]<-comp.div[up,12]
		} else {comp.div[up,15] <- "interspecific"}
}

pdf("../003.output/boxplot.divergence.by.comp.type.pdf")
with(comp.div, boxplot((genome*100)~type, ylab="percentage pairwise differences", xlab="comparison type"))
dev.off()
	
## definitely higher between species than within

#########################################################
### calculate FST?
### explaining FST by geographical distance vs. by species identity?
######################################################

pi <- read.table("../003.output/all.pop.pi.txt",stringsAsFactors=FALSE)

###########################################
### add pi to boxplots above for comparison
############################################

nt.diff <- comp.div[,c(10,15)]
pi <- as.data.frame(pi)
pi$pool.number <- rownames(pi)
pi <- merge(pi, pool.index)
pi$species <- paste(pi$species,".pi",sep="")
pi.s <- pi[,c(9,13)]
colnames(pi.s)[2]<-"type"
nt.diff <- rbind(nt.diff, pi.s)

pdf("../003.output/boxplot.divergence.and.pi.by.comp.type.pdf", width=8, height=6)
with(nt.diff, boxplot((genome*100)~type, ylab="percentage pairwise differences", xlab="comparison type"))
dev.off()

comp.fst <- rep(NA, nrow(comp.div))
for (up in 1:nrow(comp.div)){
	up.dat <- comp.div[up,]
	up.sp1 <- up.dat$Var1.pool.name
	up.sp2 <- up.dat$Var2.pool.name
	up.pi1 <- pi[pi$pool.name==up.sp1,9]
	up.pi2 <- pi[pi$pool.name==up.sp2,9]
	mean.win.pi <- mean(up.pi1, up.pi2)
	pi.between <- up.dat$genome
	fst <- (pi.between-mean.win.pi)/pi.between
	comp.fst[up]<-fst
}

comp.div$fst <- comp.fst
pdf("../003.output/boxplot.fst.by.comp.type.pdf")
with(comp.div, boxplot(fst~type, ylab="FST", xlab="comparison type"))
dev.off()


write.table(comp.div, file="../003.output/parsed.all.pairwise.divergence.txt",quote=FALSE)
save.image("pairwise.pop.divergence.Rdata")



### relative contribution of distance between populations and species status of pop combos to divergence differences?
### does this differ on chromosome four?
