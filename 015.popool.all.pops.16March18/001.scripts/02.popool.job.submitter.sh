#!/bin/bash

## automated submission of popoolation jobs
## popool.UPPOOL.sh
## DLF 04Mar15
## updated 16Mar18

cd /projects/aquilegia/004.pooled.sequencing/015.all.pool.popoolation.13March18/001.scripts

INFO_FILE=bam.poolsize.no169.txt
# tab separated file of pool bam name and pool size
SUB_FILE=01.popool.UPPOOL.sh
# need to replace UPPOOL and UPSIZE variables in this file

# changes file and qsubs
while read line
do
	name=$line
	array=(`echo $name | cut  --output-delimiter=" " -f 1-`)
	FILE=${array[0]}
	SIZE=${array[1]}
	MT=${array[2]}
	MEDIAN=${MT%?}
	TMPFILE=$FILE.tmp
	RUNFILE=popool.$FILE.sh
	echo $name
	echo $FILE
	echo $SIZE
	echo $MEDIAN
	echo $MT
	sed "s/UPPOOL/${FILE}/g" $SUB_FILE  > $TMPFILE
	mv $TMPFILE tmp2.txt
	sed "s/UPMEDIAN/${MEDIAN}/g" tmp2.txt > $TMPFILE
	sed "s/UPSIZE/${SIZE}/g" $TMPFILE > $RUNFILE
	rm $TMPFILE
	qsub $RUNFILE

done < $INFO_FILE
