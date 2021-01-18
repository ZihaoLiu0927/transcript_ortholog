#!/bin/bash
#SBATCH --mail-type=FAIL,END
#SBATCH--job-name=summaryPipeline2_MaizeWang2020.%A_%a
#SBATCH--output=/blue/mcintyre/share/transcript_distance/analysis/SLURM_LOGS/summaryPipeline2_MaizeWang2020.%A_%a.out
#SBATCH--error=/blue/mcintyre/share/transcript_distance/analysis/SLURM_LOGS/summaryPipeline2_MaizeWang2020.%A_%a.err 
#SBATCH --time=01:00:00
#SBATCH --ntasks=1
#SBATCH --mem=1gb
#SBATCH --array=1-12

module load samtools/1.10
PROJ=/blue/mcintyre/share/transcript_distance

IND_merge=$PROJ/analysis/MaizeWang2020/merge_refine_sample
OUTD_merge=$PROJ/analysis/MaizeWang2020/merge_refine_sample/summary
	mkdir -p $OUTD_merge
SAMPLE_merge=$(ls ${IND_merge}/*.merged.bam | awk -v row="$SLURM_ARRAY_TASK_ID" 'NR==row {print}')
name=$(basename ${SAMPLE_merge} ".merged.bam")

samtools view $SAMPLE_merge | wc -l > ${OUTD_merge}/${name}.counts


IND_baseline=$PROJ/analysis/MaizeWang2020/baseline_bam_subset
OUTD_baseline=$PROJ/analysis/MaizeWang2020/baseline_bam_subset/summary
        mkdir -p $OUTD_baseline
samtools view ${IND_baseline}/${name}.baseline.sam | wc -l > $OUTD_baseline/${name}.counts

IND_flair=$PROJ/analysis/MaizeWang2020/FLAIR
OUTD_flair=$PROJ/analysis/MaizeWang2020/FLAIR/summary
        mkdir -p $OUTD_flair
wc -l ${IND_flair}/${name}_flair.isoforms.bed > ${OUTD_flair}/${name}.counts









