#!/bin/sh  

cd ../../bam_files_new/

dexseq_lib_path='/bi/home/kiselevv/R/x86_64-unknown-linux-gnu-library/3.0/DEXSeq'
illumina_igenome_gtf_path='/bi/home/kiselevv/Illumina_iGenomes/Homo_sapiens/Ensembl/GRCh37/Annotation/Archives/archive-current/Genes'

python $dexseq_lib_path/python_scripts/dexseq_prepare_annotation.py $illumina_igenome_gtf_path/genes.gtf genes.gff

for i in `ls *_sn.sam`
do
	echo "python $dexseq_lib_path/python_scripts/dexseq_count.py genes.gff $i ../rna-seq-media/htseq_count_exons_new/$i.txt" | qsub -V -cwd
done

