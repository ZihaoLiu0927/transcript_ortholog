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
OUTD=${PROJ}/analysis/MaizeWang2020/lima
    mkdir -p $OUTD
LOG=${OUTD}/log
    mkdir -p $LOG
PRIMERS=$PROJ/original_data/MaizeWang2020/Wang2020_reconstructed_primers.fasta


SAMPLE=$(ls $PROJ/analysis/MaizeWang2020/ccs_merge_split/*.ccs.bam | awk -v row="$SLURM_ARRAY_TASK_ID" 'NR==row {print}')
name=$(basename $SAMPLE ".subreads.ccs.bam")
mkdir -p $OUTD/$name

lima --num-threads $SLURM_CPUS_PER_TASK \
     --isoseq \
     --peek-guess \
     --log-file ${LOG}/${name}_lima.log \
     $SAMPLE \
     $PRIMERS \
     $OUTD/$name/${name}_limademuxed.bam
