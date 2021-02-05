#!/bin/bash
#SBATCH --mail-user=zihaoliu@ufl.edu
#SBATCH --mail-type=FAIL,END
#SBATCH--job-name=sam2gtf.%A_%a
#SBATCH--output=/blue/mcintyre/share/transcript_distance/analysis/SLURM_LOGS/sam2gtf.%A_%a.out
#SBATCH--error=/blue/mcintyre/share/transcript_distance/analysis/SLURM_LOGS/sam2gtf.%A_%a.err 
#SBATCH --cpus-per-task=12
#SBATCH --time=00:30:00
#SBATCH --ntasks=1
#SBATCH --mem=2gb
#SBATCH --array=1-21


### Set Directories
PROJ=/blue/mcintyre/share/transcript_distance
DESIGN_FILE=$PROJ/design_files/design_Celegans_MaizeWang2020_noHeader.csv
DESIGN=$(cat $DESIGN_FILE | head -n $SLURM_TASK_ID | tail -n 1)
IFS=',' read -ra ARRAY <<<"$DESIGN"
SPECIES=${ARRAY[0]}
NAME=${ARRAY[2]}

OUTD=${PROJ}/analysis/${SPECIES}/sam2gtf
IND=${PROJ}/analysis/${SPECIES}/minimap2
LOG=${OUTD}/log
    mkdir -p $LOG

module load samtools/1.10
module load bedtools/2.29.2


# concerting sam to bam, then converting bam to bed12, bed12 to gtf format
# awk is used to generate a gtf format from input bed12 files
if [ $SPECIES == "Celegans" ]
then
    samtools view -b $IND/${NAME}_combinedRep.minimap2.sam > $OUTD/${NAME}.bam
elif [ $SPECIES == "MaizeWang2020" ]
then
    samtools view -b $IND/${NAME}.minimap2.sam > $OUTD/${NAME}.bam
else
    echo "invalid column in the design file!"
fi
bedtools bamtobed -split -i $OUTD/${NAME}.bam | awk -F "\t" -v source='$SPECIES' '{print $1"\t"source"\t""exon""\t"$2+1"\t"$3"\t"".""\t"$6"\t"".""\t""gene_id " "\"" $4 "\"" "; " "transcript_id " "\"" $4 "\"" ";"}' 
    > $OUTD/${NAME}.gtf
rm $OUTD/${NAME}.bam







