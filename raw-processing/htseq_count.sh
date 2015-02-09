# modules required to be activated for this script:

# python/2.7.3
# cutadapt/1.2.1

#!/bin/sh

cd ../../bam_files_new/

for i in `ls *_sn.sam`
do
	echo "htseq-count -s no -a 10 $i ~/Illumina_iGenomes/Homo_sapiens/Ensembl/GRCh37/Annotation/Archives/archive-current/Genes/genes.gtf > ../rna-seq-media/htseq_count_raw_new/$i.count" | qsub -V -cwd
done
