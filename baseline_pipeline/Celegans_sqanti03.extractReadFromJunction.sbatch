#!/bin/bash
#SBATCH --mail-user=zihaoliu@ufl.edu
#SBATCH --mail-type=FAIL,END
#SBATCH--job-name=extracReadFromJunction_Celegans.%A_%a
#SBATCH--output=/blue/mcintyre/share/transcript_distance/analysis/SLURM_LOGS/extracReadFromJunction_Celegans.%A_%a.out
#SBATCH--error=/blue/mcintyre/share/transcript_distance/analysis/SLURM_LOGS/extracReadFromJunction_Celegans.%A_%a.err
#SBATCH --time=00:10:00
#SBATCH --ntasks=1
#SBATCH --mem=3gb
#SBATCH --array=1-7

module load samtools/1.10
module load python/3.6.5

PROJ=/blue/mcintyre/share/transcript_distance/analysis/Celegans
INDLIST=$PROJ/baseline_read_filter
INDJUN=$PROJ/SQANTI3_QC_baseline
PATHSCRIPT=/blue/mcintyre/share/transcript_distance/scripts
OUTD=$PROJ/SQANTI3_filter
    mkdir -p $OUTD


translist=$(ls $INDLIST/*_read_pass_filter.txt | awk -v row="$SLURM_ARRAY_TASK_ID" 'NR==row {print}')
name=$(basename $translist "_read_pass_filter.txt")

SAMPLE_DIR=$(ls $INDLIST | awk '$0!="log" && $0!="summary"' | awk -v row="$SLURM_ARRAY_TASK_ID" 'NR==row {print}')

python ${PATHSCRIPT}/Common_extracFeatureFromList.py \
    --list $translist \
    --type class \
    --input $INDJUN/${name}/${name}_junctions.txt \
    --out $OUTD/$name.passed_junctions.txt
