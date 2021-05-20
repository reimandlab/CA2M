#!/bin/bash
for i in `seq 1 26`; do qsub -cwd -b y -P reimandlab -N shap$i -l h_rt=5000000 -l h_vmem=15g -o /dev/null -e /dev/null "source /.mounts/labs/reimandlab/private/users/oocsenas/Anaconda/conda/etc/profile.d/conda.sh; conda activate r_env; Rscript /.mounts/labs/reimandlab/private/users/oocsenas/CA2M/210317/bin/004D_run_mutsig_RFSHAP.R $i"; done