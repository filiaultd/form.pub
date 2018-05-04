#### helper functions for genome rotation scripts
#### DLF 02May18


#############################################
## proportion of sig snps overlapping a qtl
##############################################

## ss.pos is position of sig. differentiated SNPs (however you define) - a vector
## qtl.pos is start and stop position of QTL peaks - a numerical matrix
#ss.pos <- ss.n[,2]
#qtl.pos <- qtl.n[,2:3]

snp.qtl.prop <- function(ss.pos, qtl.pos){
	sq <- sapply(ss.pos, function(x){
		uq <- nrow(qtl.pos[qtl.pos[,1]<=x & qtl.pos[,2]>=x,])
		return(uq)
		})
	sq <- unlist(sq)	
	sqo <- sq[sq!=0]
	po <- length(sqo)/length(sq)
	return(po)
}

###############################################
### rotates a vector by a certain number of bp
################################################

### pos is a vector of positions to rotate
### bp.slide is the number of bases to rotate
### max.bp is the maximum positions in the genome
#pos <- ss.pos
#bp.slide <- 500000
#max.bp <- len.max

genome.rotate <- function(pos,bp.slide, max.bp){
	pos.r <- pos+bp.slide
	pos.rr <- sapply(pos.r, function(x){
		if(x>max.bp){x <- x-max.bp}
		return(x)
	})
	return(pos.rr)
}

###################################################################
## returns number of rotations in which proportion of sig snps overlapping a qtl is higher or equal to the observed proportion
#####################################################################
## ss.pos is position of sig. differentiated SNPs (however you define) - a vector
## qtl.pos is start and stop position of QTL peaks - a numerical matrix
#ss.pos <- ss.n[,2]
#qtl.pos <- qtl.n[,2:3]
## n.rotations is the number of rotations to do
## max.bp is the length in bp of the genome

qtl.snp.rotation <- function(ss.pos, qtl.pos, max.bp, n.rotations){
        ### get observed values
        obs.prop <- snp.qtl.prop(ss.pos=ss.pos, qtl.pos=qtl.pos)

        ### get rotated values
        rot.bp <- sample(1:max.bp, n.rotations, replace=FALSE)
        #print(rot.bp)
        rot.prop <- sapply(rot.bp, function(bp){
        #       print(bp)
                pos.r <- genome.rotate(pos=ss.pos, bp.slide=bp, max.bp=len.max)
                prop.r <- snp.qtl.prop(ss.pos=pos.r, qtl.pos=qtl.pos)
                return(prop.r)
                })
        rot.extreme <- rot.prop[rot.prop>=obs.prop]
        n.extreme <- length(rot.extreme)
        return(n.extreme)
        }



########################################
### get number of snps in each recombination bin
#########################################
### ss.pos is vector of positions that are sig differentiated
### rw.n is the data frame above that lists starts, stops, rec.bin of recombination windows

snp.recomb.num <- function(ss.pos, rw.n){
	sr <- sapply(ss.pos, function(x){
		#print(x)
		ur <- rw.n[rw.n$starts<=x & rw.n$stops>=x,]
		if(nrow(ur)==1){bin.out <- as.character(ur[5])
			} else {bin.out <- NA} 
		return(bin.out)
		})
	srt <- table(sr)	
	return(srt)
}





