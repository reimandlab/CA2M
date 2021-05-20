#!/bin/bash
for (( i = 1; i <= 26; i++ ))
do  
    for (( b = 1; b <= 10; b++ ))
    do  
        qsub -N boot-$i-$b -cwd -b y -P reimandlab  -l h_vmem=35g -l h_rt=300000 -o /dev/null -e /dev/null "source /.mounts/labs/reimandlab/private/users/oocsenas/Anaconda/conda/etc/profile.d/conda.sh; conda activate r_env; Rscript /.mounts/labs/reimandlab/private/users/oocsenas/CA2M/210317/bin/003C_RF_bootstrap_importances.R $i $b" 
    done 
done
