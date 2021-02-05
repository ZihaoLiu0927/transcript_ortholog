#!/bin/bash
#SBATCH --mail-type=FAIL,END
#SBATCH--job-name=CombineRep_Celegans.%A_%a
#SBATCH--output=/blue/mcintyre/share/transcript_distance/analysis/SLURM_LOGS/CombineRep_Celegans.%A_%a.out
#SBATCH--error=/blue/mcintyre/share/transcript_distance/analysis/SLURM_LOGS/CombineRep_Celegans.%A_%a.err 
#SBATCH --time=00:20:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=2gb
#SBATCH --array=1-8



PROJ=/blue/mcintyre/share/transcript_distance
IND=$PROJ/original_data/CelegansRoach2020
OUTD=$PROJ/original_data/CelegansRoach2020
DESIGN_FILE=$PROJ/design_files/design_Celegans_MaizeWang2020_noHeader.csv
DESIGN=$(cat $DESIGN_FILE | head -n $SLURM_TASK_ID | tail -n 1)
IFS=',' read -ra ARRAY <<< "$DESIGN"

SPECIES=${ARRAY[0]}
NAME=${ARRAY[2]}

cat $IND/${NAME}_rep1.fastq $IND/${NAME}_rep2.fastq > $OUTD/${NAME}_combinedRep.fastq



