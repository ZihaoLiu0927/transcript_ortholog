#!/bin/bash
#SBATCH --mail-user=zihaoliu@ufl.edu
#SBATCH --mail-type=FAIL,END
#SBATCH--job-name=Isoseq3Cluster_MaizeWang2020.%A_%a
#SBATCH--output=/blue/mcintyre/share/transcript_distance/analysis/SLURM_LOGS/Isoseq3Cluster_MaizeWang2020.%A_%a.out
#SBATCH--error=/blue/mcintyre/share/transcript_distance/analysis/SLURM_LOGS/runIsoseq3Cluster_MaizeWang2020.%A_%a.err 
#SBATCH --cpus-per-task=12
#SBATCH --time=6:00:00
#SBATCH --ntasks=1
#SBATCH --mem=16gb
#SBATCH --array=1-12

module load isoseq3/3.4.0
module load samtools/1.10

### Set Directories
PROJ=/blue/mcintyre/share/transcript_distance
OUTD=${PROJ}/analysis/MaizeWang2020/isoseq3cluster
LOG=${OUTD}/log
    mkdir -p $LOG
IND=$PROJ/analysis/MaizeWang2020/bamTocluster

SAMPLE=$(ls $IND/*.baseline.bam | awk -v row="$SLURM_ARRAY_TASK_ID" 'NR==row {print}')
name=$(basename $SAMPLE ".baseline.bam")


isoseq3 cluster \
	--num-threads $SLURM_CPUS_PER_TASK \
	--log-file $LOG/$name.log \
	--verbose \
	--use-qvs \
	$SAMPLE \
	$OUTD/${name}.isoseq3_clustered.bam
