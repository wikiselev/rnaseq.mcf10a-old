# modules required to be activated for this script:

# python/2.7.3
# cutadapt/1.2.1

#!/bin/sh  
for i in `ls *fastq.gz`
do
	echo "trim_galore -o ../trimmed $i" | qsub -V -cwd
done
