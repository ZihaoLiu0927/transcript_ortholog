#!/bin/bash
#SBATCH --mail-user=zihaoliu@ufl.edu
#SBATCH --mail-type=FAIL,END
#SBATCH--job-name=minimap2_c.elegans.%A_%a
#SBATCH--output=/blue/mcintyre/share/transcript_distance/scripts/SLURM_LOGS/minimap2_c.elegans.%A_%a.out
#SBATCH--error=/blue/mcintyre/share/transcript_distance/scripts/SLURM_LOGS/minimap2_c.elegans.%A_%a.err 
#SBATCH --cpus-per-task=12
#SBATCH --time=01:00:00
#SBATCH --ntasks=1
#SBATCH --mem=32gb
#SBATCH --array=1-14

module load minimap/2.12

### Set Directories
PROJ=/blue/mcintyre/share/transcript_distance
    mkdir -p ${PROJ}/analysis/Celegans
OUTD=${PROJ}/analysis/Celegans/minimap2
    mkdir -p $OUTD
LOG=${OUTD}/log
    mkdir -p $LOG

REF=/blue/mcintyre/share/references/Celegans_WBcel235/Caenorhabditis_elegans.WBcel235.dna.toplevel.fa

SAMPLE=$(ls ${PROJ}/data/CelegansRoach2020/*.fastq.gz | awk -v row="$SLURM_ARRAY_TASK_ID" 'NR==row {print}')
name=$(basename ${SAMPLE} ".fastq.gz")

minimap2 -t 12 -a -x splice -u f --secondary=no -C 5 $REF ${SAMPLE} \
    > $OUTD/${name}.minimap2.sam 2>$LOG/${name}.minimap2.log
