#!/bin/bash
#SBATCH --mail-user=zihaoliu@ufl.edu
#SBATCH --mail-type=FAIL,END
#SBATCH--job-name=Isoseq3Cluster_Celegans.%A_%a
#SBATCH--output=/blue/mcintyre/share/transcript_distance/analysis/SLURM_LOGS/Isoseq3Cluster_Celegans.%A_%a.out
#SBATCH--error=/blue/mcintyre/share/transcript_distance/analysis/SLURM_LOGS/Isoseq3Cluster_Celegans.%A_%a.err 
#SBATCH --cpus-per-task=12
#SBATCH --time=6:00:00
#SBATCH --ntasks=1
#SBATCH --mem=64gb
#SBATCH --array=1-8

module load isoseq3/3.4.0
module load samtools/1.10


### Set Directories
PROJ=/blue/mcintyre/share/transcript_distance
DESIGN_FILE=$PROJ/design_files/design_Celegans_MaizeWang2020_noHeader.csv
DESIGN=$(cat $DESIGN_FILE | head -n $SLURM_TASK_ID | tail -n 1)
IFS=',' read -ra ARRAY <<< "$DESIGN"
SPECIES=${ARRAY[0]}
NAME=${ARRAY[2]}

OUTD=${PROJ}/analysis/${SPECIES}/isoseq3cluster
LOG=${OUTD}/log
    mkdir -p $LOG
IND=$PROJ/analysis/${SPECIES}/baseline_bam_subset


isoseq3 cluster \
	--num-threads $SLURM_CPUS_PER_TASK \
	--log-file $LOG/$NAME.log \
	--verbose \
        --use-qvs \
	$IND/${NAME}.baseline.bam \
	$OUTD/${NAME}.isoseq3_clustered.bam
