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
OUTD=${PROJ}/analysis/MaizeWang2020/ccs
LOG=${OUTD}/log
    mkdir -p $LOG
REPORT=${OUTD}/report
    mkdir -p $REPORT


SAMPLENUM=$((($SLURM_ARRAY_TASK_ID-1)/10+1))
SAMPLE=$(ls ${PROJ}/original_data/MaizeWang2020/*.bam | awk -v row="$SAMPLENUM" 'NR==row {print}')
name=$(basename ${SAMPLE} ".subreads.bam")
chunkID=$((($SLURM_ARRAY_TASK_ID-1)%10+1))

ccs --min-passes 1 \
    --min-rq 0.9 \
    --num-threads $SLURM_CPUS_PER_TASK \
    --chunk ${chunkID}/10 \
    --log-file ${LOG}/${name}_${chunkID}.log \
    --report-file ${REPORT}/${name}_${chunkID}_report.txt \
    $SAMPLE \
    $OUTD/${name}.subreads.ccs.${chunkID}.bam

