import sys
import os

def main():
	filepath=sys.argv[1] ### whole genome sync file path
	#outfile="test.out"
	outfile=sys.argv[2]  ### output file name/path
	minCov=int(sys.argv[3])  ### minimum coverage for all pools
	#minCov=2

	if not os.path.isfile(filepath):
		print("File path {} does not exist.  Exiting".format(filepath))
		sys.exit()
	
	fileOpen = open("test.out", "w")

	poly_sites={}
	with open(filepath) as upfile:
		for lineOut in upfile:
			line=lineOut.rstrip('\n')
			#print("line {} is up next".format(line))
			line_s=line.split("\t")
			#print("splitline is {}".format(line_s))
			ref=line_s[2]
			#print("ref allele is {}".format(ref))			

			### counts of each allele
			ac = []
			tc = []
			cc = []
			gc = []
			nc = []
			delc = []
			alleles = line_s[3:13]
			#print("alleles are {}".format(alleles))
			for up_allele in alleles:
				alleleSplit(up_allele, ac, tc, cc, gc, nc, delc)			


			### sums of each allele
			asum=sum(ac)
			tsum=sum(tc)
			csum=sum(cc)
			gsum=sum(gc)
			nsum=sum(nc)
			dsum=sum(delc)
			sums=[asum, tsum, csum, gsum, nsum, dsum]
			#print(sums)						
			ssum = sum(sums)
			if ssum==0 :
				continue
			
			### sum of reference allele	
			#print("here we go!")
			bases=["A","T","C","G","N","D"]
			#baseSum=sums[bases.index(ref)]
			#print("{} sum is {}".format(ref,baseSum))

			### check sum of other alleles
			### remove if this sum is zero (i.e. position is nonpolymorphic)
			#index_remove = bases.index(ref)
			#nrsums=sums
			#del nrsums[index_remove]
			#nrSum=sum(nrsums)
			#print(nrSum)
			#if nrSum==0 :
			#	continue

			#print(line)
			### remove site that done meet min coverage threshold
			### remove nonpolymorphic sites
			### get most common allele for each pool
		
			mPos=[]
			tPos=[]
			for up_allele in alleles:
				alleleTotal(up_allele, tPos)
				alleleMax(up_allele, bases, mPos)
			

			if min(tPos)<minCov:
				continue
			
			
			#print(mPos)
			if len(set(mPos))==1:
				cPos=mPos[0]
				#print(cPos)
				nrsums=sums
                        	del nrsums[cPos]
                        	nrSum=sum(nrsums)
                        	#print(nrSum)
                        	if nrSum==0 :
                                	continue
				else:
					 fileOpen.write(lineOut)
	fileOpen.close()

				
def alleleSplit(up_allele, ac, tc, cc, gc, nc, delc):
	ua_s=up_allele.split(":")
	ua_s=map(int, ua_s)
	ac.append(ua_s[0])
	tc.append(ua_s[1])
	cc.append(ua_s[2])
	gc.append(ua_s[3])
	nc.append(ua_s[4])
	delc.append(ua_s[5])
	
def alleleMax(up_allele,bases,mPos):
	ua_s=up_allele.split(":")
	ua_s=map(int, ua_s)
	maxAl=max(ua_s)
	maxPos=ua_s.index(maxAl)
	mPos.append(maxPos)

def alleleTotal(up_allele,tPos):
	ua_s=up_allele.split(":")
        ua_s=map(int, ua_s)
        minPos=sum(ua_s)
	#print(minPos)
        tPos.append(minPos)


if __name__ == '__main__':
	main()
