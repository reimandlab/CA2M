#!/bin/bash
for (( i = 1; i <= 26; i++ ))
do  
    for (( s = 1; s <= 10; s++ ))
    do  
        qsub -N perm-$i-$s -cwd -b y -P reimandlab  -l h_vmem=35g -l h_rt=300000 -o /dev/null -e /dev/null "source /.mounts/labs/reimandlab/private/users/oocsenas/Anaconda/conda/etc/profile.d/conda.sh; conda activate r_env; Rscript /.mounts/labs/reimandlab/private/users/oocsenas/CA2M/210317/bin/004B_run_permtest_sig_importances.R $i $s" 
    done 
done
