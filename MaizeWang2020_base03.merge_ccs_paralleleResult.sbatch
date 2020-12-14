#!/bin/bash
#SBATCH --mail-user=zihaoliu@ufl.edu
#SBATCH --mail-type=FAIL,END
#SBATCH--job-name=ccsMerge_MaizeWang2020.%A_%a
#SBATCH--output=/blue/mcintyre/share/transcript_distance/analysis/SLURM_LOGS/ccsMerge_MaizeWang2020.%A_%a.out
#SBATCH--error=/blue/mcintyre/share/transcript_distance/analysis/SLURM_LOGS/ccsMerge_MaizeWang2020.%A_%a.err 
#SBATCH --time=2:00:00
#SBATCH --ntasks=1
#SBATCH --mem=3gb
#SBATCH --cpus-per-task=8
#SBATCH --array=1-15

module load samtools/1.10

### Set Directories
PROJ=/blue/mcintyre/share/transcript_distance
OUTD=${PROJ}/analysis/MaizeWang2020/ccsMerge
LOG=${OUTD}/log
    mkdir -p $LOG

SAMPLE=$(ls ${PROJ}/analysis/MaizeWang2020/ccs/*.subreads.ccs.1.bam | awk -v row="$SLURM_ARRAY_TASK_ID" 'NR==row {print}')
name=$(basename $SAMPLE ".subreads.ccs.1.bam")

samtools merge \
    --threads $SLURM_CPUS_PER_TASK \
    $OUTD/$name.subreads.ccs.bam \
    ${PROJ}/analysis/MaizeWang2020/ccs/$name.subreads.ccs.*.bam
