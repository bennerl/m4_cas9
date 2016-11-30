#!/bin/bash


##############################################################################################################
ANN=/data/bennerl/reference/
WD=/data/bennerl/m4_cas9_data/

rm ${WD}hisat2_swarm.txt
for i in ${WD}*rep*/
do
	rm ${i}*sam
	rm ${i}*txt
	FILE=${i#${WD}}
	FILE=${FILE%/}
	echo "hisat2 -x ${ANN}r6.12_dmel_noERCC --fr --dta --known-splicesite-infile ${ANN}r6.12_dmel_noERCC.splicesites.txt --novel-splicesite-outfile ${i}${FILE}_hisat_novel_splicesites.txt -p 8 -1 ${i}*1.fastq -2 ${i}*2.fastq -S ${i}${FILE}.sam" >> ${WD}hisat2_swarm.txt
done
#not --rf, rna-strandness doesnt effect alignement
swarm -g 16 -t 8 --time 01:00:00 --module hisat/2.0.4 -f ${WD}hisat2_swarm.txt
##############################################################################################################
ANN=/data/bennerl/reference/
WD=/data/bennerl/m4_cas9_data/

rm ${WD}samtools_swarm.txt
for i in ${WD}*rep*/
do
	rm ${i}*bam
	FILE=${i#$WD}
	FILE=${FILE%/}
	echo "samtools sort -@ 8 -o ${i}${FILE}.sorted.bam ${i}${FILE}.sam; samtools index ${i}${FILE}.sorted.bam; samtools sort -@ 8 -n -o ${i}${FILE}.name.sorted.bam ${i}${FILE}.sam" >> ${WD}samtools_swarm.txt
done
swarm -g 16 -t 8 --time 01:30:00 --module samtools/1.3.1 -f ${WD}samtools_swarm.txt
##############################################################################################################
ANN=/data/bennerl/reference/
WD=/data/bennerl/m4_cas9_data/
rm ${WD}stringtie_gtf_swarm.txt
for i in ${WD}*rep*/
do
	FILE=${i#${WD}}
	FILE=${FILE%/}
	echo "stringtie ${i}${FILE}.sorted.bam -o ${i}${FILE}.gtf -p 8 -G ${ANN}r6.12_dmel_noERCC.gtf" >> ${WD}stringtie_gtf_swarm.txt
done
swarm -g 16 -t 8 --time 00:30:00 --module stringtie/1.2.3 -f ${WD}stringtie_gtf_swarm.txt
##############################################################################################################
ANN=/data/bennerl/reference/
WD=/data/bennerl/m4_cas9_data/
rm ${WD}stringtie_merge_swarm.txt
rm ${WD}sample_gtf/*
for i in ${WD}*rep*/
do
	FILE=${i#${WD}}
	FILE=${FILE%/}
	cp ${i}${FILE}.gtf ${WD}sample_gtf/
done
echo "stringtie --merge -G ${ANN}r6.12_dmel_noERCC.gtf -o ${WD}sample_gtf/stringtie_merge.gtf ${WD}sample_gtf/*gtf" >> ${WD}stringtie_merge_swarm.txt
swarm -g 16 -t 8 --time 00:30:00 --module stringtie/1.2.3 -f ${WD}stringtie_merge_swarm.txt
##############################################################################################################
ANN=/data/bennerl/reference/
WD=/data/bennerl/m4_cas9_data/
rm ${WD}stringtie_ballgown_swarm.txt
for i in ${WD}*rep*/
do
	rm -r ${i}stringtie_annotated
	rm -r ${i}stringtie_merge
	FILE=${i#${WD}}
	FILE=${FILE%/}
	mkdir ${i}stringtie_annotated
	mkdir ${i}stringtie_merge
	echo "stringtie ${i}${FILE}.sorted.bam -o ${i}/stringtie_annotated/${FILE}.gtf -A ${i}/stringtie_annotated/${FILE}_gene_abundances.tab -B -e -p 8 -G ${ANN}r6.12_dmel_noERCC.gtf " >> ${WD}stringtie_ballgown_swarm.txt
	echo "stringtie ${i}${FILE}.sorted.bam -o ${i}/stringtie_merge/${FILE}.gtf -A ${i}/stringtie_merge/${FILE}_gene_abundances.tab -B -e -p 8 -G ${WD}sample_gtf/stringtie_merge.gtf" >> ${WD}stringtie_ballgown_swarm.txt
done
swarm -g 16 -t 8 --time 00:30:00 --module stringtie/1.2.3 -f ${WD}stringtie_ballgown_swarm.txt
##############################################################################################################