#!/bin/bash
#SBATCH --mail-type=FAIL,END
#SBATCH--job-name=CombineAll_MaizeWang2020.%A_%a
#SBATCH--output=/blue/mcintyre/share/transcript_distance/analysis/SLURM_LOGS/CombineAll_MaizeWang2020.%A_%a.out
#SBATCH--error=/blue/mcintyre/share/transcript_distance/analysis/SLURM_LOGS/CombineAll_MaizeWang2020.%A_%a.err 
#SBATCH --time=00:20:00
#SBATCH --ntasks=1
#SBATCH --mem=4gb


module load samtools/1.10

PROJ=/blue/mcintyre/share/transcript_distance
IND=$PROJ/analysis/MaizeWang2020/merge_refine_sample
OUTD=$PROJ/analysis/MaizeWang2020/merge_refine_sample

cat $IND/*.merged.fastq > $OUTD/MaizeWang2020.merged.fastq
