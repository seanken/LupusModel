#! /bin/bash

#$ -cwd
#$ -l h_vmem=8G
#$ -e Err/cellranger.$TASK_ID.err
#$ -o Err/cellranger.$TASK_ID.out
#$ -l h_rt=30:00:00
#$ -pe smp 8
#$ -binding linear:8
#$ -R y

#$ -t 1-9

source /broad/software/scripts/useuse
reuse Python-2.7
reuse UGER
use .cellranger-6.1.2

refDir=$1 ##directory with refernce in it

ref=${refDir}/refdata-gex-mm10-2020-A ##10X reference to use

SEEDFILE=samps.txt ##sample sheet
nam=$(awk "NR==$SGE_TASK_ID" $SEEDFILE | awk '{print $2}') ##Sample name
dir=$(awk "NR==$SGE_TASK_ID" $SEEDFILE | awk '{print $1}') ##directory with fastq in it
numCells=$(awk "NR==$SGE_TASK_ID" $SEEDFILE | awk '{print $3}') ##Expected number of cells
typ=$(awk "NR==$SGE_TASK_ID" $SEEDFILE | awk '{print $4}') ##type (3' or 5')

ID=$nam

mkdir $nam
cd $nam

cellranger count --id=samp --transcriptome=$ref --fastqs=$dir --localmem=64 --nosecondary --chemistry $typ --expect-cells $numCells --sample $nam
