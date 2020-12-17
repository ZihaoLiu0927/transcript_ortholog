#!/bin/bash
#SBATCH --mail-user=zihaoliu@ufl.edu
#SBATCH --mail-type=FAIL,END
#SBATCH--job-name=Isoseq.sam2gtf_Celegans.%A_%a
#SBATCH--output=/blue/mcintyre/share/transcript_distance/analysis/SLURM_LOGS/Isoseq.sam2gtf_Celegans.%A_%a.out
#SBATCH--error=/blue/mcintyre/share/transcript_distance/analysis/SLURM_LOGS/Isoseq.sam2gtf_Celegans.%A_%a.err 
#SBATCH --time=01:00:00
#SBATCH --ntasks=1
#SBATCH --mem=3gb
#SBATCH --array=1-7


### Set Directories
PROJ=/blue/mcintyre/share/transcript_distance
OUTD=${PROJ}/analysis/Celegans/isoseq3cluster/correctGTF
	mkdir -p $OUTD
IND=${PROJ}/analysis/Celegans/isoseq3cluster/minimap2align
LOG=${OUTD}/log
    mkdir -p $LOG

module load samtools/1.10
module load bedtools/2.29.2


SAMPLE=$(ls $IND/*.isoseq3_clustered.aligned.sam| awk -v row="$SLURM_ARRAY_TASK_ID" 'NR==row {print}')
name=$(basename $SAMPLE ".isoseq3_clustered.aligned.sam")

samtools view -b -F 2048 ${SAMPLE} > $OUTD/${name}.bam
bedtools bamtobed -split -i $OUTD/${name}.bam | awk -F "\t" '{print $1"\t""Celegans""\t""exon""\t"$2+1"\t"$3"\t"".""\t"$6"\t"".""\t""gene_id " "\"" $4 "\"" "; " "transcript_id " "\"" $4 "\"" ";"}' \
	 > $OUTD/${name}.isoseq3cluster_corrected.gtf
rm $OUTD/${name}.bam
