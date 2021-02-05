#!/bin/bash
#SBATCH --mail-user=zihaoliu@ufl.edu
#SBATCH --mail-type=FAIL,END
#SBATCH--job-name=runFLAIR_Celegans.%A_%a
#SBATCH--output=/blue/mcintyre/share/transcript_distance/analysis/SLURM_LOGS/runFLAIR_Celegans.%A_%a.out
#SBATCH--error=/blue/mcintyre/share/transcript_distance/analysis/SLURM_LOGS/runFLAIR_Celegans.%A_%a.err 
#SBATCH --time=6:00:00
#SBATCH --ntasks=1
#SBATCH --mem=8gb
#SBATCH --array=1-21

module load bedtools/2.29.2
module load flair/1.5

### Set Directories

PROJ=/blue/mcintyre/share/transcript_distance
DESIGN_FILE=$PROJ/design_files/design_Celegans_MaizeWang2020_noHeader.csv
DESIGN=$(cat $DESIGN_FILE | head -n $SLURM_TASK_ID | tail -n 1)
IFS=',' read -ra ARRAY <<< "$DESIGN"
SPECIES=${ARRAY[0]}
NAME=${ARRAY[2]}

OUTD=${PROJ}/analysis/${SPECIES}/FLAIR
BED12=$OUTD/sam2bed12
    mkdir -p $BED12
LOG=${OUTD}/log
    mkdir -p $LOG
IND=$PROJ/analysis/${SPECIES}/baseline_bam_subset
if [ $SPECIES = "Celegans" ]; 
then
    FQD=$PROJ/original_data/CelegansRoach2020
    GENOME=/blue/mcintyre/share/references/Celegans_WBcel235/Caenorhabditis_elegans.WBcel235.dna.toplevel.fa
    ANNO=/blue/mcintyre/share/references/Celegans_WBcel235/FSM_consolidation_Celegans_WBcel235/Celegans_WBcel235_FSM_consolidation_all_transcript.gtf
elif [ $SPECIES = "MaizeWang2020" ];
then
    FQD=$PROJ/analysis/MaizeWang2020/merge_refine_sample
    GENOME=/blue/mcintyre/share/references/maize_b73/b73_genome_ensemblgenomes/fasta/zea_mays/dna/Zea_mays.B73_RefGen_v4.dna.toplevel.fa
    ANNO=/blue/mcintyre/share/references/maize_b73/FSM_consolidation_maize_B73/maize_B73_FSM_consolidation_all_transcript.gtf
else
    echo "wrong source column in the design file"
fi

python $HPC_FLAIR_DIR/flair/bin/bam2Bed12.py -i $SAMPLE > ${BED12}/$NAME.bed12

python $HPC_FLAIR_DIR/flair/flair.py 23 \
    -q ${BED12}/${NAME}.bed12 \
    -g $GENOME \
    -f $ANNO \
    -r $FQD/${NAME}_combinedRep.fastq \
    -o $OUTD/${NAME}_flair \
    --keep_intermediate \
    --generate_map \
    2>$LOG/${NAME}.flair_results.log



