#!/bin/tcsh
#BSUB -J simulator[1-10]
#BSUB -n 1
#BSUB -W 01:00
#BSUB -o Output_%J_%I.out
#BSUB -e Error_%J_%I.err
#BSUB -R "rusage[mem=16GB]"

module load conda
conda activate /usr/local/usrapps/rosswhet/biostar_env

Rscript --vanilla ~/SimBreeder/cluster_submission/test.R $LSB_JOBINDEX

conda deactivate
