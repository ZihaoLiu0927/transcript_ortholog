#!/bin/bash
#SBATCH --mail-user=zihaoliu@ufl.edu
#SBATCH --mail-type=FAIL,END
#SBATCH--job-name=runIsoseq3_Celegans.%A_%a
#SBATCH--output=/blue/mcintyre/share/transcript_distance/analysis/SLURM_LOGS/runIsoseq3_Celegans.%A_%a.out
#SBATCH--error=/blue/mcintyre/share/transcript_distance/analysis/SLURM_LOGS/runIsoseq3_Celegans.%A_%a.err 
#SBATCH --cpus-per-task=12
#SBATCH --time=6:00:00
#SBATCH --ntasks=1
#SBATCH --mem=16gb
#SBATCH --array=1-7

module load isoseq3/3.4.0
module load samtools/1.10


### Set Directories
PROJ=/blue/mcintyre/share/transcript_distance
OUTD=${PROJ}/analysis/Celegans/isoseq3
    mkdir -p $OUTD
LOG=${OUTD}/log
    mkdir -p $LOG
IND=$PROJ/analysis/Celegans/minimap2

SAMPLE=$(ls $IND/*_combinedRep.minimap2.sam | awk -v row="$SLURM_ARRAY_TASK_ID" 'NR==row {print}')
name=$(basename $SAMPLE "_combinedRep.minimap2.sam")

samtools view -bS $SAMPLE > $IND/$name.bam

isoseq3 cluster \
	--num-threads $SLURM_CPUS_PER_TASK \
	--log-file $LOG/$name.log \
	--verbose \
        --use-qvs \
	$IND/$name.bam \
	$OUTD/${name}.clustered.bam
