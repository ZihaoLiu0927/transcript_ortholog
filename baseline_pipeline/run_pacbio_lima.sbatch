#!/bin/bash
#SBATCH --mail-user=zihaoliu@ufl.edu
#SBATCH --mail-type=FAIL,END
#SBATCH--job-name=lima_MaizeWang2020.%A_%a
#SBATCH--output=/blue/mcintyre/share/transcript_distance/analysis/SLURM_LOGS/lima_MaizeWang2020.%A_%a.out
#SBATCH--error=/blue/mcintyre/share/transcript_distance/analysis/SLURM_LOGS/lima_MaizeWang2020.%A_%a.err 
#SBATCH --cpus-per-task=12
#SBATCH --time=96:00:00
#SBATCH --ntasks=1
#SBATCH --mem=26gb
#SBATCH --array=1-15

module load lima/2.0.0

### Set Directories
PROJ=/blue/mcintyre/share/transcript_distance
    mkdir -p ${PROJ}/analysis/MaizeWang2020
DESIGN_FILE=$PROJ/design_files/design_MaizeWang2020_originSamples_noHeader.csv
DESIGN=$(cat $DESIGN_FILE | head -n $SLURM_ARRAY_TASK_ID | tail -n 1)
IFS=',' read -ra ARRAY <<< "$DESIGN"
NAME=${ARRAY[2]}

IND=$PROJ/analysis/MaizeWang2020/ccs_merge_split
OUTD=${PROJ}/analysis/MaizeWang2020/lima
LOG=${OUTD}/log
    mkdir -p $LOG
PRIMERS=$PROJ/original_data/MaizeWang2020/Wang2020_reconstructed_primers.fasta

mkdir -p $OUTD/$NAME

lima --num-threads $SLURM_CPUS_PER_TASK \
     --isoseq \
     --peek-guess \
     --log-file ${LOG}/${NAME}_lima.log \
     $IND/${NAME}.ccs.bam \
     $PRIMERS \
     $OUTD/$NAME/${NAME}_limademuxed.bam
