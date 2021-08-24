#! /bin/bash -l

#$ -cwd
#$ -q broad
#$ -l h_vmem=50g
#$ -e aggr3.0.err
#$ -o aggr3.0.out
#$ -l h_rt=18:00:00
#$ -l os=RedHat7

reuse Python-2.7
reuse UGER
reuse .cellranger-3.0.2
cellranger aggr --id=EXPR --csv=/path/to/aggr_sampleSheet.csv --normalize=mapped