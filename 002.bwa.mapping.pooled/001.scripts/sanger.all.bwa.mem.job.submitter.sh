#!/bin/bash

## automated submission of bwa mem mapping jobs
## sanger.all.bwa.mem.UPSPECIES.sh
## DLF 07 July 14
## for new genome: 24 March 15

for FILE in 302GBAAXX_1_8281_20091109_seq.interleaved.fastq  302GBAAXX_7_8289_20091109_seq.interleaved.fastq  8283_sequence.interleaved.fastq  8434_sequence.interleaved.fastq  302GBAAXX_2_8282_20091109_seq.interleaved.fastq  302GBAAXX_8_8290_20091109_seq.interleaved.fastq  8284_sequence.interleaved.fastq  8435_sequence.interleaved.fastq  302GBAAXX_3_8283_20091109_seq.interleaved.fastq  81KHFABXX_6_31-08-2011.interleaved.fastq         8285_sequence.interleaved.fastq  C0HWGACXX_1_20120713B_20120724.interleaved.fastq  302GBAAXX_4_8284_20091109_seq.interleaved.fastq  81KHFABXX_7_31-08-2011.interleaved.fastq  8288_sequence.interleaved.fastq  s_1_sequence.interleaved.fastq  302GBAAXX_5_8285_20091109_seq.interleaved.fastq  8281_sequence.interleaved.fastq  8289_sequence.interleaved.fastq  s_2_sequence.interleaved.fastq  302GBAAXX_6_8288_20091109_seq.interleaved.fastq  8282_sequence.interleaved.fastq  8290_sequence.interleaved.fastq  s_6_sequence.interleaved.fastq

do sed "s/UPSPECIES/${FILE}/g" sanger.all.bwa.mem.UPSPECIES.sh  > sanger.all.bwa.mem.$FILE.sh.tmp
IDUP=`echo $FILE | cut -d"." -f1`
NAMEUP=${IDUP:0:4}
sed "s/JOBID/$IDUP/g" sanger.all.bwa.mem.$FILE.sh.tmp > sanger.all.bwa.mem.$FILE.sh.tmp2
sed "s/NAMEID/$NAMEUP/g" sanger.all.bwa.mem.$FILE.sh.tmp2 > sanger.all.bwa.mem.$FILE.sh
qsub sanger.all.bwa.mem.$FILE.sh
rm sanger.all.bwa.mem.$FILE.sh.tmp
rm sanger.all.bwa.mem.$FILE.sh.tmp2
done


 
