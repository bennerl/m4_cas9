#!/bin/bash
ANN=/data/bennerl/reference/
WD=/data/bennerl/m4_cas9_data/

rm ${WD}hisat2_swarm.txt
for i in ${WD}*/
do
	FILE=${i#$WD}
	FILE=${FILE%/}
	echo "hisat2 -x ${ANN}r6.12_dmel_noERCC --fr --dta --known-splicesite-infile ${ANN}r6.12_dmel_noERCC.splicesites.txt --novel-splicesite-outfile ${i}${FILE}_hisat_novel_splicesites.txt -p 8 -1 ${i}*1.fastq -2 ${i}*2.fastq -S ${i}${FILE}.sam" >> ${WD}/hisat2_swarm.txt
done
swarm -g 16 -t 8 --time 01:00:00 --module hisat/2.0.4 -f ${WD}hisat2_swarm.txt




