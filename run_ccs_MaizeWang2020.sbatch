#!/bin/bash
#SBATCH --mail-user=zihaoliu@ufl.edu
#SBATCH --mail-type=FAIL,END
#SBATCH--job-name=ccs_MaizeWang2020.%A_%a
#SBATCH--output=/blue/mcintyre/share/transcript_distance/scripts/SLURM_LOGS/ccs_MaizeWang2020.%A_%a.out
#SBATCH--error=/blue/mcintyre/share/transcript_distance/scripts/SLURM_LOGS/ccs_MaizeWang2020.%A_%a.err 
#SBATCH --time=96:00:00
#SBATCH --ntasks=1
#SBATCH --mem=16gb
#SBATCH --cpus-per-task=12
#SBATCH --array=1-15

### isoseq3 analysis - library processing
### have mixed sample with 3 replicates


module load pbccs/5.0.0 

### Set Directories
PROJ=/blue/mcintyre/share/transcript_distance
    mkdir -p ${PROJ}/analysis/MaizeWang2020
OUTD=${PROJ}/analysis/MaizeWang2020/ccs
LOG=${OUTD}/log
    mkdir -p $LOG
REPORT=${OUTD}/report
    mkdir -p $REPORT


SAMPLE=$(ls ${PROJ}/data/MaizeWang2020/*.bam | awk -v row="$SLURM_ARRAY_TASK_ID" 'NR==row {print}')
name=$(basename ${SAMPLE} ".subreads.bam")
 ccs --min-passes 1 --num-threads 12 \
     --log-file ${LOG}/${name}.log \
     --report-file ${REPORT}/${name}_report.txt \
     $SAMPLE \
     $OUTD/${name}.subreads.ccs.bam

