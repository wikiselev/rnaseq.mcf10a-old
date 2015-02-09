#!/bin/sh  
for i in `ls`
do
	# echo $i
	#statements
	echo "fastqc $i -o fastqc_out" | qsub -V -cwd -N $i
done

rm *fastq.gz.*
rm *.e*
rm *.o*