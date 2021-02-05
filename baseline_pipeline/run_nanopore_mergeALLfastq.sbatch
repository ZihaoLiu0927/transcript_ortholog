#!/bin/bash
#SBATCH --mail-type=FAIL,END
#SBATCH--job-name=CombineAll_Celegans.%A_%a
#SBATCH--output=/blue/mcintyre/share/transcript_distance/analysis/SLURM_LOGS/CombineAll_Celegans.%A_%a.out
#SBATCH--error=/blue/mcintyre/share/transcript_distance/analysis/SLURM_LOGS/CombineAll_Celegans.%A_%a.err 
#SBATCH --time=00:20:00
#SBATCH --ntasks=1
#SBATCH --mem=4gb


module load samtools/1.10

PROJ=/blue/mcintyre/share/transcript_distance
IND=$PROJ/original_data/CelegansRoach2020
OUTD=$PROJ/original_data/CelegansRoach2020

cat $IND/*_rep*.fastq > $OUTD/Celegans_combinedRep.fastq
