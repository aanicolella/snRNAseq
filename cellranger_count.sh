#! /bin/bash -l

#$ -cwd
#$ -l h_vmem=8G
#$ -e /path/to/mkfastq/REPLICATE/cellranger-3.0.err
#$ -o /path/to/mkfastq/REPLICATE/cellranger-3.0.out
#$ -l h_rt=48:00:00
#$ -pe smp 8
#$ -binding linear:8
#$ -R y
#$ -l os=RedHat7

reuse Python-2.7
reuse UGER
reuse .cellranger-3.0.2
cellranger count --id=REPLICATE --transcriptome=/path/to/reference/transcriptome --expect-cells=NUMCELLS --fastqs=/path/to/mkfastq/REPLICATE --sample=REPLICATE --chemistry=SC3Pv3 --localmem=64
