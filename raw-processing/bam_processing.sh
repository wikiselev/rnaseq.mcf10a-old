# modules required to be activated for this script:

# python/2.7.3
# cutadapt/1.2.1

#!/bin/sh  
for i in `ls ../../mapped_new/`
do
	echo "sh samtools.sh $i" | qsub -V -cwd
done
