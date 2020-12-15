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
IND=$PROJ/analysis/MaizeWang2020/isoseqRefine
OUTD=${PROJ}/analysis/MaizeWang2020/mergedSample
    mkdir -p $OUTD/sortInputFile 
LOG=${OUTD}/log
    mkdir -p $LOG
GENOTYPE=$PROJ/data/MaizeWang2020/Wang2020_reconstructed_primers.fasta


module load samtools/1.10
module load bedtools/2.29.2

for i in `ls $IND/*/*.refined.bam`
do
    name=$(basename $i ".bam") 
    samtools sort $i > $OUTD/sortInputFile/$name.sorted.bam
done


for i in `grep -Po "(?<=>).*(?=_3p)" $GENOTYPE`
do
    echo $i >> $LOG/file_merge.log
    echo "merged file:" >> $LOG/file_merge.log
    ls $OUTD/sortInputFile/*.$i.refined.sorted.bam | awk '{print}' >> $LOG/file_merge.log
    samtools merge $OUTD/$i.merged.bam $OUTD/sortInputFile/*.${i}.refined.sorted.bam
done

#for i in `ls $OUTD/*.merged.bam`
#do
#    name=$(basename $i ".merged.bam")
#    bedtools bamtofastq -i $i \
#	                -fq $OUTD/$name.merged.fastq
#done

rm -rf $OUTD/sortInputFile
#rm $OUTD/*.merged.bam



