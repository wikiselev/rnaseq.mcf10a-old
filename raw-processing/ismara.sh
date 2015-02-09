for i in `ls -d *_trimmed`
do
    echo "cp $i/accepted_hits.bam ../ismara/$i.bam" | qsub -V -cwd
done

