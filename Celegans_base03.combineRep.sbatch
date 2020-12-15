#!/bin/bash
#SBATCH --mail-type=FAIL,END
#SBATCH--job-name=CombineRep_Celegans.%A_%a
#SBATCH--output=/blue/mcintyre/share/transcript_distance/analysis/SLURM_LOGS/CombineRep_Celegans.%A_%a.out
#SBATCH--error=/blue/mcintyre/share/transcript_distance/analysis/SLURM_LOGS/CombineRep_Celegans.%A_%a.err 
#SBATCH --time=00:20:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=3gb
#SBATCH --array=1-7


module load samtools/1.10

PROJ=/blue/mcintyre/share/transcript_distance
IND=$PROJ/analysis/Celegans/minimap2
OUTD=$PROJ/analysis/Celegans/minimap2

SAMPLE=$(ls ${IND}/*_rep1.minimap2.sam | awk -v row="$SLURM_ARRAY_TASK_ID" 'NR==row {print}')
name=$(basename $SAMPLE "_rep1.minimap2.sam")

grep -v "@" $IND/${name}_rep2.minimap2.sam > $OUTD/${name}.temp.sam
cat $SAMPLE $OUTD/${name}.temp.sam | samtools sort -@ $SLURM_CPUS_PER_TASK -o $OUTD/${name}_combinedRep.minimap2.sam

rm $OUTD/${name}.temp.sam


