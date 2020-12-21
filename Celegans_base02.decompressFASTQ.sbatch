#!/bin/bash
#SBATCH --mail-user=zihaoliu@ufl.edu
#SBATCH --mail-type=FAIL,END
#SBATCH--job-name=decompressFASTQ_Celegans.%A_%a
#SBATCH--output=/blue/mcintyre/share/transcript_distance/analysis/SLURM_LOGS/decompressFASTQ_Celegans.%A_%a.out
#SBATCH--error=/blue/mcintyre/share/transcript_distance/analysis/SLURM_LOGS/decompressFASTQ_Celegans.%A_%a.err 
#SBATCH --time=1:00:00
#SBATCH --ntasks=1
#SBATCH --mem=4gb
#SBATCH --array=1-14

PROJ=/blue/mcintyre/share/transcript_distance/original_data/CelegansRoach2020

SAMPLE=$(ls $PROJ/*.fastq.gz | awk -v row="$SLURM_ARRAY_TASK_ID" 'NR==row {print}')
name=$(basename $SAMPLE ".gz")

gzip -dc $SAMPLE > $PROJ/$name
