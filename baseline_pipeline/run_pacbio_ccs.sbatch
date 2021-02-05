#!/bin/bash
#SBATCH --mail-user=zihaoliu@ufl.edu
#SBATCH --mail-type=FAIL,END
#SBATCH--job-name=ccs_MaizeWang2020.%A_%a
#SBATCH--output=/blue/mcintyre/share/transcript_distance/analysis/SLURM_LOGS/ccs_MaizeWang2020.%A_%a.out
#SBATCH--error=/blue/mcintyre/share/transcript_distance/analysis/SLURM_LOGS/ccs_MaizeWang2020.%A_%a.err 
#SBATCH --time=96:00:00
#SBATCH --ntasks=1
#SBATCH --mem=16gb
#SBATCH --cpus-per-task=12
#SBATCH --array=1-21

### isoseq3 analysis - library processing
### have mixed sample with 3 replicates


module load pbccs/5.0.0 

### Set Directories
PROJ=/blue/mcintyre/share/transcript_distance
    mkdir -p ${PROJ}/analysis/MaizeWang2020
DESIGN_FILE=$PROJ/design_files/design_MaizeWang2020_originSamples_noHeader.csv
DESIGN=$(cat $DESIGN_FILE | head -n $SLURM_CPUS_PER_TASK | tail -n 1)
IFS=',' read -ra ARRAY <<< "$DESIGN"
NAME=${ARRAY[2]}

IND=$PROJ/original_data/MaizeWang2020
OUTD=${PROJ}/analysis/MaizeWang2020/ccs
LOG=${OUTD}/log
    mkdir -p $LOG
REPORT=${OUTD}/report
    mkdir -p $REPORT

ccs --min-passes 1 --num-threads $SLURM_CPUS_PER_TASK \
     --log-file ${LOG}/${NAME}.log \
     --report-file ${REPORT}/${NAME}_report.txt \
     $IND/${NAME}.subreads.bam \
     $OUTD/${NAME}.subreads.ccs.bam

