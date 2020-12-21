#!/bin/bash
#SBATCH --mail-type=FAIL,END
#SBATCH--job-name=CombineRep_Celegans.%A_%a
#SBATCH--output=/blue/mcintyre/share/transcript_distance/analysis/SLURM_LOGS/CombineRep_Celegans.%A_%a.out
#SBATCH--error=/blue/mcintyre/share/transcript_distance/analysis/SLURM_LOGS/CombineRep_Celegans.%A_%a.err 
#SBATCH --time=00:20:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=2gb
#SBATCH --array=1-7


module load samtools/1.10

PROJ=/blue/mcintyre/share/transcript_distance
IND=$PROJ/original_data/CelegansRoach2020
OUTD=$PROJ/original_data/CelegansRoach2020

SAMPLE=$(ls ${IND}/*_rep1.fastq | awk -v row="$SLURM_ARRAY_TASK_ID" 'NR==row {print}')
name=$(basename $SAMPLE "_rep1.fastq")

cat $SAMPLE $IND/${name}_rep2.fastq > $OUTD/${name}_combinedRep.fastq



