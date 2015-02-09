#!/bin/sh
for i in `ls ../../mapped_new_hg19/`
do
	echo "cp ../../mapped_new_hg19/$i/accepted_hits.bam ../../bam_files_new_hg19/$i.bam" | qsub -V -cwd
done
