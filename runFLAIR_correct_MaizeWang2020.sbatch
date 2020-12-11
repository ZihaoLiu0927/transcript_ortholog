#!/bin/bash
#SBATCH --mail-user=zihaoliu@ufl.edu
#SBATCH --mail-type=FAIL,END
#SBATCH--job-name=runFLAIR_MaizeWang2020.%A_%a
#SBATCH--output=/blue/mcintyre/share/transcript_distance/analysis/SLURM_LOGS/runFLAIR_MaizeWang2020.%A_%a.out
#SBATCH--error=/blue/mcintyre/share/transcript_distance/analysis/SLURM_LOGS/runFLAIR_MaizeWang2020.%A_%a.err 
#SBATCH --cpus-per-task=12
#SBATCH --time=96:00:00
#SBATCH --ntasks=1
#SBATCH --mem=8gb
#SBATCH --array=1-12

module load bedtools/2.29.2
#require kerneltree installed

### Set Directories
PROJ=/blue/mcintyre/share/transcript_distance
OUTD=${PROJ}/analysis/MaizeWang2020/FLAIR/FLAIRcorrect
    mkdir -p $OUTD
BED12=$OUTD/sam2bed12
    mkdir -p $BED12
LOG=${OUTD}/log
    mkdir -p $LOG
IND=$PROJ/analysis/MaizeWang2020/minimap2
SCRIPT=$PROJ/tools/flair
GENOME=/blue/mcintyre/share/references/maize_b73/b73_genome_ensemblgenomes/fasta/zea_mays/dna/Zea_mays.B73_RefGen_v4.dna/toplevel.fa.gz

SAMPLE=$(ls $IND | awk -v row="$SLURM_ARRAY_TASK_ID" 'NR==row {print}')
name=$(basename $SAMPLE ".minimap2.sam")

${SCRIPT}/bin/bam2Bed12.py -i $SAMPLE > ${BED12}/$name.bed12

python ${SCRIPT}/flair.py correct \
-q ${BED12}/$name.bed \
-g $GENOME \
-o $OUTD/flair



