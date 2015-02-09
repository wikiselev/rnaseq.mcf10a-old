#!/bin/sh  
# sort by name, convert to SAM for htseq-count
samtools sort -n ../../mapped_new/$1/accepted_hits.bam ../../bam_files_new/$1_sn
samtools view -h -o ../../bam_files_new/$1_sn.sam ../../bam_files_new/$1_sn.bam

# sort by position and index for IGV
# samtools sort ../../mapped_new/$1/accepted_hits.bam ../../bam_files_new/$1_sp
# samtools view -h -o ../../bam_files_new/$1_sp.sam ../../bam_files_new/$1_sp.bam
# samtools index ../../bam_files_new/$1_sp.sam
