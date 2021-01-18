#!/bin/bash
#SBATCH --mail-user=zihaoliu@ufl.edu
#SBATCH --mail-type=FAIL,END
#SBATCH--job-name=bam2fastq_MaizeWang2020.%A_%a
#SBATCH--output=/blue/mcintyre/share/transcript_distance/analysis/SLURM_LOGS/bam2fastq_MaizeWang2020.%A_%a.out
#SBATCH--error=/blue/mcintyre/share/transcript_distance/analysis/SLURM_LOGS/bam2fastq_MaizeWang2020.%A_%a.err 
#SBATCH --time=5:00:00
#SBATCH --ntasks=1
#SBATCH --mem=4gb
#SBATCH --array=1-15

module load bedtools/2.29.2
#require kerneltree installed

### Set Directories
PROJ=/blue/mcintyre/share/transcript_distance
OUTD=${PROJ}/data/MaizeWang2020
    mkdir -p $OUTD
LOG=${OUTD}/log
    mkdir -p $LOG

SAMPLE=$(ls $OUTD/*.subreads.bam | awk -v row="$SLURM_ARRAY_TASK_ID" 'NR==row {print}')
name=$(basename $SAMPLE ".subreads.bam")

bedtools bamtofastq \
	-i $SAMPLE \
	-fq $OUTD/$name.fastq
