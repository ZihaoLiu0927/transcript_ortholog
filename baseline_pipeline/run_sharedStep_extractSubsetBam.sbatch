#!/bin/bash
#SBATCH --mail-user=zihaoliu@ufl.edu
#SBATCH --mail-type=FAIL,END
#SBATCH--job-name=extracTransFromAlign_Celegans.%A_%a
#SBATCH--output=/blue/mcintyre/share/transcript_distance/analysis/SLURM_LOGS/extracTransFromAlign_Celegans.%A_%a.out
#SBATCH--error=/blue/mcintyre/share/transcript_distance/analysis/SLURM_LOGS/extracTransFromAlign_Celegans.%A_%a.err
#SBATCH --time=01:10:00
#SBATCH --ntasks=1
#SBATCH --mem=32gb
#SBATCH --array=1-21

module load samtools/1.10

PROJ=/blue/mcintyre/share/transcript_distance
DESIGN_FILE=$PROJ/design_files/design_Celegans_MaizeWang2020_noHeader.csv
DESIGN=$(cat $DESIGN_FILE | head -n $SLURM_TASK_ID | tail -n 1)
IFS=',' read -ra ARRAY <<< "$DESIGN"
SPECIES=${ARRAY[0]}
NAME=${ARRAY[2]}

INDLIST=$PROJ/analysis/${SPECIES}/baseline_read_filter
INDSAM=$PROJ/analysis/${SPECIES}/minimap2
if [ $SPECIES=="Celegans" ];
then
    INPUT=$INDSAM/${NAME}_combinedRep.minimap2.sam
elif
then
    INPUT=$INDSAM/${NAME}.minimap2.sam
fi
PATHSCRIPT=/blue/mcintyre/share/transcript_distance/scripts
OUTD=$PROJ/analysis/${SPECIES}/baseline_bam_subset
    mkdir -p $OUTD

module load python/3.6.5
python ${PATHSCRIPT}/Common_extracFeatureFromList.py \
    --list $INDLIST/${NAME}_read_pass_filter.txt \
    --type sam \
    --input $INPUT \
    --out $OUTD/$NAME.baseline.sam

samtools view -bS $OUTD/$NAME.baseline.sam > $OUTD/$NAME.baseline.bam
