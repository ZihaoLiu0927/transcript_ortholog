#!/bin/bash
#SBATCH --mail-user=zihaoliu@ufl.edu
#SBATCH --mail-type=FAIL,END
#SBATCH--job-name=extracTransFromAlign_MaizeWang2020.%A_%a
#SBATCH--output=/blue/mcintyre/share/transcript_distance/analysis/SLURM_LOGS/extracTransFromAlign_MaizeWang2020.%A_%a.out
#SBATCH--error=/blue/mcintyre/share/transcript_distance/analysis/SLURM_LOGS/extracTransFromAlign_MaizeWang2020.%A_%a.err 
#SBATCH --time=00:10:00
#SBATCH --ntasks=1
#SBATCH --mem=8gb
#SBATCH --array=1-12

module load samtools/1.10


PROJ=/blue/mcintyre/share/transcript_distance/analysis/MaizeWang2020
INDLIST=$PROJ/baseline_read_filter
INDSAM=$PROJ/minimap2
PATHSCRIPT=/blue/mcintyre/share/transcript_distance/scripts
OUTD=$PROJ/bamTocluster
    mkdir -p $OUTD

translist=$(ls $INDLIST/*.merged_read_pass_filter.txt | awk -v row="$SLURM_ARRAY_TASK_ID" 'NR==row {print}')
name=$(basename $translist ".merged_read_pass_filter.txt")

module load python/3.6.5
python ${PATHSCRIPT}/Common_extracTransFromAlign.py --list $translist --sam $INDSAM/$name.merged.fastq.minimap2.sam --out $OUTD/$name.baseline.sam

samtools view -bS $OUTD/$name.baseline.sam > $OUTD/$name.baseline.bam
