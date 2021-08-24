#! /bin/bash

#$ -cwd
#$ -q broad
#$ -l h_vmem=20g
#$ -e /path/to/data/ambientRNA/REP/REP.err
#$ -o /path/to/data/ambientRNA/REP/REP.out
#$ -l h_rt=120:00:00
#$ -l os=RedHat7
#$ -pe smp 16
#$ -binding linear:16
#$ -R y

source /broad/software/scripts/useuse

use .java-jdk-1.8.0_181-x86-64
use .hdf5-1.8.16
use .anaconda3-5.3.1

echo Init!

echo Time to test!!

source activate cellbender_v2
dir=/path/to/cellRanger/count/outputs/REP/outs/raw_feature_bc_matrix
output=/path/to/data/ambientRNA/REP/REP.results.h5
cellbender remove-background --input $dir --output $output --expected-cells NUMCELLS --epochs 300 --total-droplets-included MAXCELLS --training-fraction .2 --z-dim 200 --z-layers 1000
conda deactivate
