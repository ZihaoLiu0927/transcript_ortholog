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


### Set Directories
PROJ=/blue/mcintyre/share/transcript_distance
    mkdir -p ${PROJ}/analysis/MaizeWang2020
OUTD=${PROJ}/analysis/MaizeWang2020/FLAIR/FLAIRcollapse
    mkdir -p $OUTD
LOG=${OUTD}/log
    mkdir -p $LOG
IND=$PROJ/analysis/MaizeWang2020/FLAIR/FLAIRcorrect
FA=
SCRIPT=$PROJ/script/FLAIR/bin
GENOME=/blue/mcintyre/share/references/maize_b73/b73_genome_ensemblgenomes/fasta/zea_mays/dna/Zea_mays.B73_RefGen_v4.dna/toplevel.fa.gz


SAMPLE=$(ls $IND/*.bed | awk -v row="$SLURM_ARRAY_TASK_ID" 'NR==row {print}')
name=$(basename $SAMPLE ".bed")

python flair.py collapse \
-t $SLURM_CPUS_PER_TASK \
-g $GENOME \
-r $FA \
-q $IND/$SAMPLE \
-o $OUD/flair.collapse





