#!/bin/bash
#SBATCH --mail-user=zihaoliu@ufl.edu
#SBATCH --mail-type=FAIL,END
#SBATCH--job-name=Isoseq.sam2gtf.%A_%a
#SBATCH--output=/blue/mcintyre/share/transcript_distance/analysis/SLURM_LOGS/Isoseq.sam2gtf.%A_%a.out
#SBATCH--error=/blue/mcintyre/share/transcript_distance/analysis/SLURM_LOGS/Isoseq.sam2gtf.%A_%a.err 
#SBATCH --time=01:00:00
#SBATCH --ntasks=1
#SBATCH --mem=6gb
#SBATCH --array=1-21


### Set Directories
PROJ=/blue/mcintyre/share/transcript_distance
DESIGN_FILE=$PROJ/design_files/design_Celegans_MaizeWang2020_noHeader.csv
DESIGN=$(cat $DESIGN_FILE | head -n $SLURM_TASK_ID | tail -n 1)
IFS=',' read -ra ARRAY <<< "$DESIGN"
SPECIES=${ARRAY[0]}
NAME=${ARRAY[2]}

OUTD=${PROJ}/analysis/${SPECIES}/isoseq3cluster/isoseqClusterGTF
IND=${PROJ}/analysis/${SPECIES}/isoseq3cluster/minimap2align
LOG=${OUTD}/log
    mkdir -p $LOG

module load samtools/1.10
module load bedtools/2.29.2

# concerting sam to bam, then converting bam to bed12, bed12 to gtf format
# awk is used to generate a gtf format from input bed12 files
samtools view -b -F 2048 $IND/${NAME}.isoseq3_clustered.aligned.sam > $OUTD/${NAME}.bam
bedtools bamtobed -split -i $OUTD/${NAME}.bam | awk -F "\t" '{print $1"\t""Celegans""\t""exon""\t"$2+1"\t"$3"\t"".""\t"$6"\t"".""\t""gene_id " "\"" $4 "\"" "; " "transcript_id " "\"" $4 "\"" ";"}' \
	 > $OUTD/${NAME}.isoseq3cluster_corrected.gtf
rm $OUTD/${NAME}.bam
