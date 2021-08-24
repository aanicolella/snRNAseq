#! /bin/bash -l

#$ -cwd
#$ -q broad
#$ -l h_vmem=50g
#$ -e mkfastq3.0.err
#$ -o mkfastq3.0.out
#$ -l h_rt=18:00:00
#$ -l os=RedHat7

reuse .bcl2fastq2-v2.20.0
reuse Python-2.7
reuse UGER
reuse .cellranger-3.0.2
cellranger mkfastq --run=/path/to/rawData --samplesheet=/path/to/sampleSheet.csv 
 --output-dir=/path/to/output/directory/mkfastq/ --localmem=96 --jobmode=local