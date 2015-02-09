# modules required to be activated for this script:

# python/2.7.3
# cutadapt/1.2.1

#!/bin/sh  
for i in `ls ../../trimmed/*fq.gz | xargs -n 1 basename`
do
	arr=(${i//./ })
	mkdir ../../mapped_new_hg19/$arr
	echo "tophat2 -o ../../mapped_new_hg19/$arr ~/Illumina_iGenomes/Homo_sapiens/UCSC/hg19/Sequence/Bowtie2Index/genome ../../trimmed/$i" | qsub -V -cwd
done
