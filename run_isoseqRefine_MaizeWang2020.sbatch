#!/bin/bash
#SBATCH --mail-user=zihaoliu@ufl.edu
#SBATCH --mail-type=FAIL,END
#SBATCH--job-name=iso-refine_MaizeWang2020.%A_%a
#SBATCH--output=/blue/mcintyre/share/transcript_distance/analysis/SLURM_LOGS/iso-refine_MaizeWang2020.%A_%a.out
#SBATCH--error=/blue/mcintyre/share/transcript_distance/analysis/SLURM_LOGS/iso-refine_MaizeWang2020.%A_%a.err 
#SBATCH --cpus-per-task=24
#SBATCH --time=96:00:00
#SBATCH --ntasks=1
#SBATCH --mem=10gb
#SBATCH --array=1-15

module load isoseq3/3.4.0

### Set Directories
PROJ=/blue/mcintyre/share/transcript_distance
OUTD=${PROJ}/analysis/MaizeWang2020/isoseqRefine
LOG=${OUTD}/log
    mkdir -p $LOG
PRIMERS=$PROJ/data/MaizeWang2020/Wang2020_reconstructed_primers.fasta


SAMPLE_DIR=$(ls ${PROJ}/analysis/MaizeWang2020/lima | awk '$0!="log" && $0!="summary"' | awk -v row="$SLURM_ARRAY_TASK_ID" 'NR==row {print}')
mkdir -p $OUTD/$SAMPLE_DIR



for i in `ls $PROJ/analysis/MaizeWang2020/lima/$SAMPLE_DIR/*.bam`
do
	genotype=$(basename $i "_3p.bam" | grep -Po "(?<=--).*")
	isoseq3 refine \
		--num-threads $SLURM_CPUS_PER_TASK \
		--require-polya \
		$i \
	        $PRIMERS \
		$OUTD/$SAMPLE_DIR/${SAMPLE_DIR}.${genotype}.refined.bam
done		

