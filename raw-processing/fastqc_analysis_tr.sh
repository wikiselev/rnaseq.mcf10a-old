#!/bin/sh  
for i in `ls *.fq.gz`
do
	# echo $i
	#statements
	echo "fastqc $i -o fastqc_out" | qsub -V -cwd -N $i
done

rm *fq.gz.*
