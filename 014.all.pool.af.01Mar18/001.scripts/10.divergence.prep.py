import sys
import os

#### script to parse an rsync input to generate prep files for divergence/pi calculations
#### output file one: chr pos 0/1(pass coverage threshold for each pool)
#### output file two: only biallelic sites (sync file)



def main():
	filepath=sys.argv[1] ### whole genome sync file path
	#outfile="test.out"
	outfile=sys.argv[2]  ### prefix for output file name/path
	minCov=int(sys.argv[3])  ### minimum coverage for pools
	#minCov=2
	minPass=int(sys.argv[4]) ### minimum number of pools that pass min cov threshold

	if not os.path.isfile(filepath):
		print("File path {} does not exist.  Exiting".format(filepath))
		sys.exit()
	
	biFILE = open((outfile + ".biallelic.out"), "w")
	goodFILE = open((outfile + ".goodsites.out"),"w")

	with open(filepath) as upfile:
		for lineOut in upfile:
			line=lineOut.rstrip('\n')
			line_s=line.split("\t")
			ref=line_s[2]

			### counts of each allele across all pools
			ac = []
			tc = []
			cc = []
			gc = []
			nc = []
			delc = []
			alleles = line_s[3:13]
			for up_allele in alleles:
				alleleSplit(up_allele, ac, tc, cc, gc, nc, delc)			


			### sums of each allele
			asum=sum(ac)
			tsum=sum(tc)
			csum=sum(cc)
			gsum=sum(gc)
			nsum=sum(nc)
			dsum=sum(delc)
			sums=[asum, tsum, csum, gsum]
			#print(sums)				
			
			### remove sites with deletions in any pool
			if dsum > 0 :
				continue
		
			### remove sites with no coverage in all pools
			ssum = sum(sums)
			if ssum==0 :
				continue

			### keep sites that are non-poly or bialleic only
			bpSum = sums[0:4]
			bpNZ = [num for num in bpSum if num > 0]
			if len(bpNZ) > 2:
				continue
						
			# need to assess if min number pools pass coverage threshold
			### cov of each pool
                        alleles = line_s[3:13]
                        pCov=[]
                        for up_allele in alleles:
                                poolCov(up_allele, pCov)

			passCov = []
			for upCov in pCov:
				if upCov > minCov:
					passCov.append(1)
				else:
					passCov.append(0)
			
			### remove sites where not min cov in min number of pools
			pcSum = sum(passCov)
			if pcSum < minPass:
				continue

			
			## output these good sites in compact format
			passCov = map(str, passCov)	
			passCovCat = ''.join(passCov)

                        ### generate output line
                        pos=line_s[0:2]
                        pos.append(passCovCat)

                        ## write to output file
                        for element in pos:
                                goodFILE.write(element)
                                goodFILE.write("\t")
                       	goodFILE.write("\n")

			

			### output good bialleic sites only
                        if len(bpNZ) == 2:
                                biFILE.write(lineOut)


	biFILE.close()
	goodFILE.close()
				
def alleleSplit(up_allele, ac, tc, cc, gc, nc, delc):
	ua_s=up_allele.split(":")
	ua_s=map(int, ua_s)
	ac.append(ua_s[0])
	tc.append(ua_s[1])
	cc.append(ua_s[2])
	gc.append(ua_s[3])
	nc.append(ua_s[4])
	delc.append(ua_s[5])
	

def poolCov(up_allele,pCov):
        ua_s=up_allele.split(":")
        ua_s=map(int, ua_s)
        ua_sum=sum(ua_s)
        pCov.append(ua_sum)



#def alleleMax(up_allele,bases,mPos):
#	ua_s=up_allele.split(":")
#	ua_s=map(int, ua_s)
#	maxAl=max(ua_s)
#	maxPos=ua_s.index(maxAl)
#	mPos.append(maxPos)

#def alleleTotal(up_allele,tPos):
#	ua_s=up_allele.split(":")
#        ua_s=map(int, ua_s)
#        minPos=sum(ua_s)
#	print(minPos)
#        tPos.append(minPos)


if __name__ == '__main__':
	main()
