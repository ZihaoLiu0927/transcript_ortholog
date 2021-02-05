#!/bin/bash
#SBATCH --mail-user=zihaoliu@ufl.edu
#SBATCH --mail-type=FAIL,END
#SBATCH--job-name=decompressFASTQ_Celegans.%A_%a
#SBATCH--output=/blue/mcintyre/share/transcript_distance/analysis/SLURM_LOGS/decompressFASTQ_Celegans.%A_%a.out
#SBATCH--error=/blue/mcintyre/share/transcript_distance/analysis/SLURM_LOGS/decompressFASTQ_Celegans.%A_%a.err 
#SBATCH --time=1:00:00
#SBATCH --ntasks=1
#SBATCH --mem=4gb
#SBATCH --array=1-7

PROJ=/blue/mcintyre/share/transcript_distance
DESIGN_FILE=$PROJ/design_files/design_Celegans_MaizeWang2020_noHeader.csv
data=$PROJ/original_data/CelegansRoach2020

DESIGN=$(cat $DESIGN_FILE | head -n $SLURM_TASK_ID | tail -n 1)
IFS=',' read -ra ARRAY <<< "$DESIGN"

SPECIES=${ARRAY[0]}
NAME=${ARRAY[2]}

gzip -dc $data/${NAME}_rep1.fastq.gz > $data/${NAME}_rep1.fastq
gzip -dc $data/${NAME}_rep2.fastq.gz > $data/${NAME}_rep2.fastq


