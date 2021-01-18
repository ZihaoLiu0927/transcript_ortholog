#!/bin/bash
#SBATCH --mail-user=zihaoliu@ufl.edu
#SBATCH --mail-type=FAIL,END
#SBATCH--job-name=runFLAIR_MaizeWang2020.%A_%a
#SBATCH--output=/blue/mcintyre/share/transcript_distance/analysis/SLURM_LOGS/runFLAIR_MaizeWang2020.%A_%a.out
#SBATCH--error=/blue/mcintyre/share/transcript_distance/analysis/SLURM_LOGS/runFLAIR_MaizeWang2020.%A_%a.err 
#SBATCH --time=1:00:00
#SBATCH --ntasks=1
#SBATCH --mem=4gb
#SBATCH --array=1-12

module load bedtools/2.29.2
module load flair/1.5

### Set Directories
PROJ=/blue/mcintyre/share/transcript_distance
OUTD=${PROJ}/analysis/MaizeWang2020/FLAIR
BED12=$OUTD/sam2bed12
    mkdir -p $BED12
LOG=${OUTD}/log
    mkdir -p $LOG
IND=$PROJ/analysis/MaizeWang2020/baseline_bam_subset
FQD=$PROJ/analysis/MaizeWang2020/merge_refine_sample
GENOME=/blue/mcintyre/share/references/maize_b73/b73_genome_ensemblgenomes/fasta/zea_mays/dna/Zea_mays.B73_RefGen_v4.dna.toplevel.fa
ANNO=/blue/mcintyre/share/references/maize_b73/FSM_consolidation_maize_B73/maize_B73_FSM_consolidation_all_transcript.gtf



SAMPLE=$(ls $IND/*.baseline.sam | awk -v row="$SLURM_ARRAY_TASK_ID" 'NR==row {print}')
name=$(basename $SAMPLE ".baseline.sam")

python $HPC_FLAIR_DIR/flair/bin/bam2Bed12.py -i $SAMPLE > ${BED12}/$name.bed12

python $HPC_FLAIR_DIR/flair/flair.py 23 \
    -q ${BED12}/${name}.bed12 \
    -g $GENOME \
    -f $ANNO \
    -r $FQD/${name}.merged.fastq \
    -o $OUTD/${name}_flair \
    --keep_intermediate \
    --generate_map \
    2>$LOG/${name}.flair_results.log



