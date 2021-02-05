#!/bin/bash
#SBATCH --mail-user=zihaoliu@ufl.edu
#SBATCH --mail-type=FAIL,END
#SBATCH--job-name=minimap2.%A_%a
#SBATCH--output=/blue/mcintyre/share/transcript_distance/analysis/SLURM_LOGS/minimap2.%A_%a.out
#SBATCH--error=/blue/mcintyre/share/transcript_distance/analysis/SLURM_LOGS/minimap2.%A_%a.err 
#SBATCH --cpus-per-task=12
#SBATCH --time=01:00:00
#SBATCH --ntasks=1
#SBATCH --mem=8gb
#SBATCH --array=1-21

module load minimap/2.12
module load samtools/1.10

### Set Directories
PROJ=/blue/mcintyre/share/transcript_distance
    mkdir -p ${PROJ}/analysis/${SPECIES}
DESIGN_FILE=$PROJ/design_files/design_Celegans_MaizeWang2020_noHeader.csv
DESIGN=$(cat $DESIGN_FILE | head -n $SLURM_TASK_ID | tail -n 1)
IFS=',' read -ra ARRAY <<< "$DESIGN"
SPECIES=${ARRAY[0]}
NAME=${ARRAY[2]}
OUTD=${PROJ}/analysis/${SPECIES}/minimap2
    mkdir -p $OUTD
LOG=${OUTD}/log
    mkdir -p $LOG
if [ $SPECIES = "Celegans" ]
then
    IND=$PROJ/original_data/CelegansRoach2020
    REF=/blue/mcintyre/share/references/Celegans_WBcel235/Caenorhabditis_elegans.WBcel235.dna.toplevel.fa
elif [ $SPECIES = "MaizeWang2020" ]
then
    IND=${PROJ}/analysis/MaizeWang2020/merge_refine_sample
    REF=/blue/mcintyre/share/references/maize_b73/b73_genome_ensemblgenomes/fasta/zea_mays/dna/Zea_mays.B73_RefGen_v4.dna.toplevel.fa.gz
else
    echo "invalid source column in the design file!"
fi


minimap2 -t $SLURM_CPUS_PER_TASK -a -x splice -u f --secondary=no -C 5 $REF $IND/${NAME}_combinedRep.fastq \
    > $OUTD/${NAME}.unfiltered.minimap2.sam 2>$LOG/${NAME}.minimap2.log

samtools view -h -F 2052 $OUTD/${NAME}.unfiltered.minimap2.sam > $OUTD/${NAME}.minimap2.sam
rm $OUTD/${NAME}.unfiltered.minimap2.sam
