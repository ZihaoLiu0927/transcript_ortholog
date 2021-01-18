#!/bin/bash
#SBATCH --mail-user=zihaoliu@ufl.edu
#SBATCH --mail-type=FAIL,END
#SBATCH--job-name=minimap2_maize2020.%A_%a
#SBATCH--output=/blue/mcintyre/share/transcript_distance/analysis/SLURM_LOGS/minimap2_maize2020.%A_%a.out
#SBATCH--error=/blue/mcintyre/share/transcript_distance/analysis/SLURM_LOGS/minimap2_maize2020.%A_%a.err 
#SBATCH --cpus-per-task=12
#SBATCH --time=01:00:00
#SBATCH --ntasks=1
#SBATCH --mem=32gb
#SBATCH --array=1-12

module load minimap/2.12
module load samtools/1.10


### Set Directories
PROJ=/blue/mcintyre/share/transcript_distance
OUTD=${PROJ}/analysis/MaizeWang2020/minimap2
LOG=${OUTD}/log
    mkdir -p $LOG
REF=/blue/mcintyre/share/references/maize_b73/b73_genome_ensemblgenomes/fasta/zea_mays/dna/Zea_mays.B73_RefGen_v4.dna.toplevel.fa.gz


SAMPLE=$(ls ${PROJ}/analysis/MaizeWang2020/merge_refine_sample/*.fastq | awk -v row="$SLURM_ARRAY_TASK_ID" 'NR==row {print}')
name=$(basename ${SAMPLE} ".merged.fastq")

minimap2 -t $SLURM_CPUS_PER_TASK -a -x splice -u f --secondary=no -C 5 $REF ${SAMPLE} \
    > $OUTD/${name}.unfiltered.minimap2.sam 2>$LOG/${name}.unfiltered.minimap2.log

samtools view -h -F 2052 $OUTD/${name}.unfiltered.minimap2.sam > $OUTD/${name}.minimap2.sam
rm $OUTD/${name}.unfiltered.minimap2.sam
