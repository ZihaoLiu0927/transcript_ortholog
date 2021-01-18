#!/bin/bash
#SBATCH --mail-type=FAIL,END
#SBATCH--job-name=summaryPipeline_MaizeWang2020.%A_%a
#SBATCH--output=/blue/mcintyre/share/transcript_distance/analysis/SLURM_LOGS/summaryPipeline_MaizeWang2020.%A_%a.out
#SBATCH--error=/blue/mcintyre/share/transcript_distance/analysis/SLURM_LOGS/summaryPipeline_MaizeWang2020.%A_%a.err 
#SBATCH --time=01:00:00
#SBATCH --ntasks=1
#SBATCH --mem=1gb
#SBATCH --array=1-15

module load samtools/1.10
PROJ=/blue/mcintyre/share/transcript_distance


OUTD_origin=$PROJ/original_data/MaizeWang2020/summary
	mkdir -p $OUTD_origin
SAMPLE_origin=$(ls ${PROJ}/original_data/MaizeWang2020/*.subreads.bam | awk -v row="$SLURM_ARRAY_TASK_ID" 'NR==row {print}')
echo $SAMPLE_origin " sample"
echo $name " name"
name=$(basename ${SAMPLE_origin} ".subreads.bam")
samtools view $SAMPLE_origin | wc -l > ${OUTD_origin}/${name}.counts



OUTD_ccs=$PROJ/analysis/MaizeWang2020/ccs_merge_split/summary
	mkdir -p $OUTD_ccs
SAMPLE_ccs=${PROJ}/analysis/MaizeWang2020/ccs_merge_split/${name}.subreads.ccs.bam
samtools view $SAMPLE_ccs | wc -l > ${OUTD_ccs}/$name.counts



OUTD_refine=$PROJ/analysis/MaizeWang2020/isoseqRefine/summary
	mkdir -p $OUTD_refine
for i in {1..3}
do
	SAMPLE_refine=$(ls ${PROJ}/analysis/MaizeWang2020/isoseqRefine/$name/${name}.*.refined.bam | awk -v row="$i" 'NR==row {print}')
	name_geno=$(basename ${SAMPLE_refine} ".refined.bam")
	samtools view $SAMPLE_refine | wc -l > ${OUTD_refine}/${name_geno}.counts
done





