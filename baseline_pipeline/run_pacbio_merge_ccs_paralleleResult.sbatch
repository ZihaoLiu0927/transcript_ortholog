#!/bin/bash
#SBATCH --mail-user=zihaoliu@ufl.edu
#SBATCH --mail-type=FAIL,END
#SBATCH--job-name=ccsMerge_MaizeWang2020.%A_%a
#SBATCH--output=/blue/mcintyre/share/transcript_distance/analysis/SLURM_LOGS/ccsMerge_MaizeWang2020.%A_%a.out
#SBATCH--error=/blue/mcintyre/share/transcript_distance/analysis/SLURM_LOGS/ccsMerge_MaizeWang2020.%A_%a.err 
#SBATCH --time=2:00:00
#SBATCH --ntasks=1
#SBATCH --mem=6gb
#SBATCH --cpus-per-task=8
#SBATCH --array=1-15

module load samtools/1.10

### Set Directories
PROJ=/blue/mcintyre/share/transcript_distance
DESIGN_FILE=$PROJ/design_files/design_MaizeWang2020_originSamples_noHeader.csv
DESIGN=$(cat $DESIGN_FILE | head -n $SLURM_ARRAY_TASK_ID | tail -n 1)
IFS=',' read -ra ARRAY <<< "$DESIGN"
NAME=${ARRAY[2]}

IND=${PROJ}/analysis/MaizeWang2020/ccs
OUTD=${PROJ}/analysis/MaizeWang2020/ccsMerge
LOG=${OUTD}/log
    mkdir -p $LOG

samtools merge \
    --threads $SLURM_CPUS_PER_TASK \
    $OUTD/$NAME.subreads.ccs.bam \
    $IND/${NAME}.subreads.ccs.*.bam
