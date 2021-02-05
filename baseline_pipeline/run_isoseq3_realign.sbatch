#!/bin/bash
#SBATCH --mail-user=zihaoliu@ufl.edu
#SBATCH --mail-type=FAIL,END
#SBATCH--job-name=isoseq3_realign_Celegans.%A_%a
#SBATCH--output=/blue/mcintyre/share/transcript_distance/analysis/SLURM_LOGS/isoseq3.minimap2_realign_Celegans.%A_%a.out
#SBATCH--error=/blue/mcintyre/share/transcript_distance/analysis/SLURM_LOGS/isoseq3.minimap2_realign_Celegans.%A_%a.err 
#SBATCH --cpus-per-task=12
#SBATCH --time=1:00:00
#SBATCH --ntasks=1
#SBATCH --mem=16gb
#SBATCH --array=2

module load minimap/2.12
module load bedtools/2.29.2

### Set Directories
PROJ=/blue/mcintyre/share/transcript_distance
DESIGN_FILE=$PROJ/design_files/design_Celegans_MaizeWang2020_noHeader.csv
DESIGN=$(cat $DESIGN_FILE | head -n $SLURM_TASK_ID | tail -n 1)
IFS=',' read -ra ARRAY <<< "$DESIGN"
SPECIES=${ARRAY[0]}
NAME=${ARRAY[2]}

IND=${PROJ}/analysis/${SPECIES}/isoseq3cluster
OUTD=${PROJ}/analysis/${SPECIES}/isoseq3cluster/minimap2align
LOG=${OUTD}/log
    mkdir -p $LOG
if [ $SPECIES = "Celegans" ]
then
    REF=/blue/mcintyre/share/references/Celegans_WBcel235/Caenorhabditis_elegans.WBcel235.dna.toplevel.fa.gz
elif [ $SPECIES = "MaizeWang2020" ]
then
    REF=/blue/mcintyre/share/references/maize_b73/b73_genome_ensemblgenomes/fasta/zea_mays/dna/Zea_mays.B73_RefGen_v4.dna.toplevel.fa.gz
else
    echo "invalid source column in the design file!"
fi

minimap2 -t $SLURM_CPUS_PER_TASK -a -x splice -u f --secondary=no -C 5 $REF $IND/${NAME}.isoseq3_clustered.hq.fasta.gz \
    > $OUTD/${NAME}.aligned.sam 2>$LOG/${NAME}.aligned.log
