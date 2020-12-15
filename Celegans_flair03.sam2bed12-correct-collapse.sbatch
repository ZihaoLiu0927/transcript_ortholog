#!/bin/bash
#SBATCH --mail-user=zihaoliu@ufl.edu
#SBATCH --mail-type=FAIL,END
#SBATCH--job-name=runFLAIR_Celegans.%A_%a
#SBATCH--output=/blue/mcintyre/share/transcript_distance/analysis/SLURM_LOGS/runFLAIR_Celegans.%A_%a.out
#SBATCH--error=/blue/mcintyre/share/transcript_distance/analysis/SLURM_LOGS/runFLAIR_Celegans.%A_%a.err 
#SBATCH --time=1:00:00
#SBATCH --ntasks=1
#SBATCH --mem=4gb
#SBATCH --array=1-7

module load bedtools/2.29.2
module load flair/1.5

### Set Directories
PROJ=/blue/mcintyre/share/transcript_distance
OUTD=${PROJ}/analysis/Celegans/FLAIR
    mkdir -p $OUTD
BED12=$OUTD/sam2bed12
    mkdir -p $BED12
LOG=${OUTD}/log
    mkdir -p $LOG
IND=$PROJ/analysis/Celegans/minimap2
FQD=$PROJ/data/CelegansRoach2020
GENOME=/blue/mcintyre/share/references/Celegans_WBcel235/Caenorhabditis_elegans.WBcel235.dna.toplevel.fa
ANNO=/blue/mcintyre/share/references/Celegans_WBcel235/FSM_consolidation_Celegans_WBcel235/Celegans_WBcel235_FSM_consolidation_all_transcript.gtf



SAMPLE=$(ls $IND/*_combinedRep.minimap2.sam | awk -v row="$SLURM_ARRAY_TASK_ID" 'NR==row {print}')
name=$(basename $SAMPLE "_combinedRep.minimap2.sam")

python $HPC_FLAIR_DIR/flair/bin/bam2Bed12.py -i $SAMPLE > ${BED12}/$name.bed12

python $HPC_FLAIR_DIR/flair/flair.py 23 \
    -q ${BED12}/${name}.bed12 \
    -g $GENOME \
    -f $ANNO \
    -r $FQD/${name}_rep1.fastq \
       $FQD/${name}_rep2.fastq \
    -o $OUTD/${name}_flair \
    --keep_intermediate \
    2>$LOG/${name}.flair_results.log



