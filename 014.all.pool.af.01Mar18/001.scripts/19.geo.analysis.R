### geographical relationships between pools
### effect of geography vs "species" effect
### DLF 26April18

setwd("/Volumes/aquilegia/004.pooled.sequencing/014.all.pool.af.01Mar18/002.input/")

#p.loc <- read.csv("pool.sample.summary.table.csv",stringsAsFactors=FALSE) 
p.loc <- read.csv("sample.summary.table.30April18.csv", stringsAsFactors=FALSE)  ### made sure pool locations and short names match Christos' paper and Scott's email
p.loc <- p.loc[p.loc$df.pool.number%in%c(1,6,9)==FALSE,]
p.loc$pool.name <- gsub(" ","",p.loc$pool.name)

div <- read.table("../003.output/parsed.all.pairwise.divergence.txt")

pool.index <- read.table("pool.order.index.txt",stringsAsFactors=FALSE, header=TRUE)
pool.index$pool.name <- gsub("2","Trail",pool.index$pool.name)  ## now pool names match up between p.loc and pool.index
pool.index$pool.name <- gsub("RockCreek/","",pool.index$pool.name)  ## now pool names match up between p.loc and pool.index
pool.index <- pool.index[,-7]
pool.index <- merge(pool.index, p.loc, by=c("pool.name"))
pool.index <- pool.index[order(pool.index$order),]

write.table(pool.index, "combined.pool.index.30April18.txt", quote=FALSE)
####################################################################################
#### first make a map of pool locations to use in interpretation
#####################################################################################

library(maps)
library(maptools)
data(us.cities)
library(rgdal)

pools.f <- p.loc[p.loc$species=="f",]
pools.p <- p.loc[p.loc$species=="p",]

postscript("../003.output/pool.map.ps",width=8,height=8)
par(xpd=NA)
map("state","california",fill=FALSE,new=TRUE,lwd=1,xlim=c(min(p.loc$long)-0.6,max(p.loc$long)+0.4),ylim=c(min(p.loc$lat)-0.3,max(p.loc$lat)+0.1))
map.cities(us.cities, country="CA",minpop=100000,label=TRUE,pch=19,cex=2.2)
points(pools.f$long,pools.f$lat,pch=19,col="red",cex=2)
points(pools.p$long,pools.p$lat,pch=19,col="blue2",cex=2)
map.scale(x=min(p.loc$long)-0.6,ratio=FALSE,relwidth=0.2,cex=1.3)
legend(-118.7,36.6,pch=19,col=c("red","blue2"),legend=c("A. formosa","A. pubescens"),cex=1.4, text.font=3)
to.left <- p.loc[-c(1,3,4,7),]
text(to.left$short.name, x=to.left$long, y=to.left$lat, pos=2)
to.right <- p.loc[c(1,3,4,7),]
text(to.right$short.name, x=to.right$long, y=to.right$lat, pos=4)

dev.off()

### this is pretty boring map, but enough for me to use as a reference



####################################################################################
#### next get geographic distance between population pairs
#####################################################################################library(geosphere)
library(geosphere)
pool.index <- pool.index[order(pool.index$order),]

#vector of two numbers, a matrix of 2 columns (first one is longitude, second is latitude) or a SpatialPoints* object
#up.dist <- (distHaversine(p.loc[1,2:3],p.loc[2,2:3]))/1000

comp.dist <- rep(NA,nrow(div))
for(up in 1:nrow(div)){
	pop1 <- div[up,1]
	pop2 <- div[up,2]
	up.dist <- ((distHaversine(pool.index[pop1,7:8],pool.index[pop2,7:8])))/1000
	comp.dist[up] <- up.dist
}

hist(comp.dist)

div <- cbind(div, comp.dist)

with(div, boxplot(comp.dist~type,ylab="km"))  ### ok, kinda makes sense...

lm1 <- lm(genome~type,data=div)  ###1.483e-14 ***
lm2 <- lm(genome~comp.dist, data=div)  ### 0.6318
lm3 <- lm(genome~comp.dist, data=div[div$type=="pubescens",]) ###0.5877
lm4 <- lm(genome~comp.dist, data=div[div$type=="formosa",])  ### 0.887
lm5 <- lm(genome~comp.dist, data=div[div$type=="interspecific",])  ### 0.8626

lm6 <- lm(fst~type,data=div)  ###2.891e-11***
lm7 <- lm(Chr_04~type,data=div)  ###2.921e-06 ***


### should also try this with "genetic" distance - calculated in previous script

#save.image("19.geo.analysis.Rdata")

######################################################################
#### mantel tests with snp distance
######################################################

load("snp.dist.Rdata")

### make distances into a matrix
phys.dist <- matrix(0, ncol=10, nrow=10)
for (up in 1:nrow(div)){
	up.dat <- div[up,]
	phys.dist[up.dat[,1], up.dat[,2]] <- up.dat[,17]
	phys.dist[up.dat[,2], up.dat[,1]] <- up.dat[,17]

}
phys.dist <- as.dist(phys.dist)
library(ade4)
mt1 <- mantel.rtest(phys.dist, snp.dist, nrepet=9999)  #Simulated p-value: 0.914
### Not significant


### do this by type matrix (0=within, 1=between)
type.dist <- matrix(0, ncol=10, nrow=10)
for (up in 1:nrow(div)){
	up.dat <- div[up,]
	if(up.dat$type=="interspecific"){
		type.dist[up.dat[,1], up.dat[,2]] <- 1
		type.dist[up.dat[,2], up.dat[,1]] <- 1
	}
}
type.dist <- as.dist(type.dist)
mt2 <- mantel.rtest(type.dist, snp.dist, nrepet=9999)  #Simulated p-value: 0.0075
## yes significant


###################
### what about just within each individual species?
###########################
phys.dist <- matrix(0, ncol=10, nrow=10)
for (up in 1:nrow(div)){
	up.dat <- div[up,]
	phys.dist[up.dat[,1], up.dat[,2]] <- up.dat[,17]
	phys.dist[up.dat[,2], up.dat[,1]] <- up.dat[,17]

}
snp.dist <- as.matrix(snp.dist)

ppd <- phys.dist[pool.index$species.x=="pubescens",pool.index$species.x=="pubescens"]
psd <- snp.dist[pool.index$species.x=="pubescens",pool.index$species.x=="pubescens"]
fpd <- phys.dist[pool.index$species.x=="formosa",pool.index$species.x=="formosa"]
fsd <- snp.dist[pool.index$species.x=="formosa",pool.index$species.x=="formosa"]

phys.dist <- as.dist(phys.dist)
snp.dist <- as.dist(snp.dist)
psd <- as.dist(psd)
fsd <- as.dist(fsd)
ppd <- as.dist(ppd)
fpd <- as.dist(fpd)


pman <- mantel.rtest(psd, ppd, nrepet=9999) # pval= 0.6171 
fman <- mantel.rtest(fsd, fpd, nrepet=9999) # pval=0.9782

###################################
### mantel with log distance
#####################################

phys.dist.log <- log(phys.dist)
mt3 <- mantel.rtest(phys.dist.log, snp.dist, nrepet=9999)  #Simulated p-value: 0.9697

#########################################
### using FST/1-FST
#########################################
fst.dist <- matrix(0, ncol=10, nrow=10)
for (up in 1:nrow(div)){
	up.dat <- div[up,]
	fst.dist[up.dat[,1], up.dat[,2]] <- (up.dat[,16])/(1-up.dat[,16])
	fst.dist[up.dat[,2], up.dat[,1]] <- (up.dat[,16])/(1-up.dat[,16])
}

fst.dist <- as.dist(fst.dist)
mt4 <- mantel.rtest(phys.dist, fst.dist, nrepet=9999)  #Simulated p-value: 0.8614
mt5 <- mantel.rtest(phys.dist.log, fst.dist, nrepet=9999)  #Simulated p-value: 0.9041


##############################################################
### partial tests with both genetic distance and species
#############################################################
#library(ncf)
#pm1 <- partial.mantel.test(M1 = phys.dist, M2 = type.dist, M3 = snp.dist, resamp = 500)
#doesn't work, but doesn't make much sense anyway... the species differences are what matters

############################ 
### plot of div by distance
############################
#quartz(width=6, height=6)
pdf("../003.output/divergence.by.distance.pdf",width=6, height=6)
with(div, plot(genome*100~comp.dist, ylab="percent pairwise differences", xlab="pool separation (km)", type="n"))
with(div[div$type=="interspecific",], points(genome*100~comp.dist, col="green",pch=19))
with(div[div$type=="pubescens",], points(genome*100~comp.dist, col="blue",pch=19))
with(div[div$type=="formosa",], points(genome*100~comp.dist, col="red",pch=19))
legend(x="bottomright", col=c("green", "red","blue"), legend=c("interspecific","formosa", "pubescens"), pch=19, pt.cex=1.4)
dev.off()


save.image("19.geo.analysis.Rdata")





