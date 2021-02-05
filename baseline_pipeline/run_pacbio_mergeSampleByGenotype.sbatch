#!/bin/bash
#SBATCH --mail-user=zihaoliu@ufl.edu
#SBATCH --mail-type=FAIL,END
#SBATCH--job-name=mergelimaResults_MaizeWang2020.%A_%a
#SBATCH--output=/blue/mcintyre/share/transcript_distance/analysis/SLURM_LOGS/mergelimaResults_MaizeWang2020.%A_%a.out
#SBATCH--error=/blue/mcintyre/share/transcript_distance/analysis/SLURM_LOGS/mergelimaResults_MaizeWang2020.%A_%a.err 
#SBATCH --cpus-per-task=12
#SBATCH --time=01:00:00
#SBATCH --ntasks=1
#SBATCH --mem=3gb


### Set Directories
PROJ=/blue/mcintyre/share/transcript_distance
DESIGN_FILE=$PROJ/design_files/design_MaizeWang2020_originSamples_noHeader.csv
SAMPLES=$(cat $DESIGN_FILE | cut -d , -f 3)

IND=$PROJ/analysis/MaizeWang2020/isoseqRefine
OUTD=${PROJ}/analysis/MaizeWang2020/merge_refine_sample
    mkdir -p $OUTD/sortInputFile 
LOG=${OUTD}/log
    mkdir -p $LOG
GENOTYPE=$PROJ/original_data/MaizeWang2020/Wang2020_reconstructed_primers.fasta

DESIGN_FILE_AFTER=$PROJ/design_files/design_Celegans_MaizeWang2020_noHeader.csv
GENOTYPES=$(cat $DESIGN_FILE_AFTER | awk -F , '$1=="MaizeWang2020" {print}' | cut -d , -f 3)

module load samtools/1.10
module load bedtools/2.29.2

## sort each input file to be combined
for i in $SAMPLES
do
    for j in $(ls $IND/$i/${i}.*.bam)
    do
	name=$(basename $j ".bam") 
	samtools sort $j > $OUTD/sortInputFile/$name.sorted.bam
    done
done

##  
for geno in $GENOTYPE
do
    echo $geno >> $LOG/file_merge.log
    echo "merged file:" >> $LOG/file_merge.log
    ls $OUTD/sortInputFile/*.${geno}.refined.sorted.bam | awk '{print}' >> $LOG/file_merge.log
    samtools merge $OUTD/${geno}.bam $OUTD/sortInputFile/*.${geno}.refined.sorted.bam
    bedtools bamtofastq -i $OUTD/${geno}.bam \
	                -fq $OUTD/${geno}.fastq
done

rm -rf $OUTD/sortInputFile
rm $OUTD/*.merged.bam



