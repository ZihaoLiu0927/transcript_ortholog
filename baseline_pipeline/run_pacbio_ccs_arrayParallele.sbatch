#!/bin/bash
#SBATCH --mail-user=zihaoliu@ufl.edu
#SBATCH --mail-type=FAIL,END
#SBATCH--job-name=ccs_MaizeWang2020.%A_%a
#SBATCH--output=/blue/mcintyre/share/transcript_distance/analysis/SLURM_LOGS/ccs_MaizeWang2020.%A_%a.out
#SBATCH--error=/blue/mcintyre/share/transcript_distance/analysis/SLURM_LOGS/ccs_MaizeWang2020.%A_%a.err 
#SBATCH --time=96:00:00
#SBATCH --ntasks=1
#SBATCH --mem=4gb
#SBATCH --cpus-per-task=10
#SBATCH --array=1-150

### isoseq3 analysis - library processing
### have mixed sample with 3 replicates


module load pbccs/5.0.0 

### Set Directories
PROJ=/blue/mcintyre/share/transcript_distance
    mkdir -p ${PROJ}/analysis/MaizeWang2020
SAMPLENUM=$((($SLURM_ARRAY_TASK_ID-1)/10+1))
chunkID=$((($SLURM_ARRAY_TASK_ID-1)%10+1))

DESIGN_FILE=$PROJ/design_files/design_MaizeWang2020_originSamples_noHeader.csv
DESIGN=$(cat $DESIGN_FILE | head -n $SAMPLENUM | tail -n 1)
IFS=',' read -ra ARRAY <<< "$DESIGN"
NAME=${ARRAY[2]}

IND=$PROJ/original_data/MaizeWang2020
OUTD=${PROJ}/analysis/MaizeWang2020/ccs
LOG=${OUTD}/log
    mkdir -p $LOG
REPORT=${OUTD}/report
    mkdir -p $REPORT

ccs --min-passes 1 \
    --min-rq 0.9 \
    --num-threads $SLURM_CPUS_PER_TASK \
    --chunk ${chunkID}/10 \
    --log-file ${LOG}/${NAME}_${chunkID}.log \
    --report-file ${REPORT}/${NAME}_${chunkID}_report.txt \
    $IND/${NAME}.subreads.bam \
    $OUTD/${NAME}.subreads.ccs.${chunkID}.bam

