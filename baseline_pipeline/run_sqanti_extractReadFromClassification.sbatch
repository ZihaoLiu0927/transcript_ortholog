#!/bin/bash
#SBATCH --mail-user=zihaoliu@ufl.edu
#SBATCH --mail-type=FAIL,END
#SBATCH--job-name=extracReadFromClassification.%A_%a
#SBATCH--output=/blue/mcintyre/share/transcript_distance/analysis/SLURM_LOGS/extracReadFromClassification.%A_%a.out
#SBATCH--error=/blue/mcintyre/share/transcript_distance/analysis/SLURM_LOGS/extracReadFromClassification.%A_%a.err
#SBATCH --time=00:10:00
#SBATCH --ntasks=1
#SBATCH --mem=3gb
#SBATCH --array=1-21

module load samtools/1.10
module load python/3.6.5

PROJ=/blue/mcintyre/share/transcript_distance
DESIGN_FILE=$PROJ/design_files/design_Celegans_MaizeWang2020_noHeader.csv
DESIGN=$(cat $DESIGN_FILE | head -n $SLURM_TASK_ID | tail -n 1)
IFS=',' read -ra ARRAY <<< "$DESIGN"
SPECIES=${ARRAY[0]}
NAME=${ARRAY[2]}

INDLIST=$PROJ/analysis/$SPECIES/baseline_read_filter
INDCLS=$PROJ/analysis/$SPECIES/SQANTI3_QC_baseline
PATHSCRIPT=$PROJ/scripts
OUTD=$PROJ/SQANTI3_filter
    mkdir -p $OUTD

python ${PATHSCRIPT}/Common_extracFeatureFromList.py \
    --list $INDLUST/${NAME}_read_pass_filter.txt \
    --type class \
    --input $INDCLS/${NAME}/${NAME}_classification.txt \
    --out $OUTD/$NAME.passed_classification.txt
