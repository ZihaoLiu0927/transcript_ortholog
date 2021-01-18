#!/bin/bash
#SBATCH --mail-user=zihaoliu@ufl.edu
#SBATCH --mail-type=FAIL,END
#SBATCH--job-name=isoseq3_realign_maize2020.%A_%a
#SBATCH--output=/blue/mcintyre/share/transcript_distance/analysis/SLURM_LOGS/isoseq3.minimap2_realign_maize2020.%A_%a.out
#SBATCH--error=/blue/mcintyre/share/transcript_distance/analysis/SLURM_LOGS/isoseq3.minimap2_realign_maize2020.%A_%a.err 
#SBATCH --cpus-per-task=12
#SBATCH --time=01:00:00
#SBATCH --ntasks=1
#SBATCH --mem=16gb
#SBATCH --array=1-12

module load minimap/2.12
module load bedtools/2.29.2

### Set Directories
PROJ=/blue/mcintyre/share/transcript_distance
IND=${PROJ}/analysis/MaizeWang2020/isoseq3cluster
OUTD=${PROJ}/analysis/MaizeWang2020/isoseq3cluster/minimap2align
LOG=${OUTD}/log
    mkdir -p $LOG
REF=/blue/mcintyre/share/references/maize_b73/b73_genome_ensemblgenomes/fasta/zea_mays/dna/Zea_mays.B73_RefGen_v4.dna.toplevel.fa.gz

#convert bam file to fastq file

SAMPLE=$(ls $IND/*.isoseq3_clustered.hq.fasta.gz | awk -v row="$SLURM_ARRAY_TASK_ID" 'NR==row {print}')
name=$(basename ${SAMPLE} ".hq.fasta.gz")


minimap2 -t $SLURM_CPUS_PER_TASK -a -x splice -u f --secondary=no -C 5 $REF $SAMPLE \
    > $OUTD/${name}.aligned.sam 2>$LOG/${name}.aligned.log
