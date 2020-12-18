#!/bin/bash
#SBATCH --mail-user=zihaoliu@ufl.edu
#SBATCH --mail-type=FAIL,END
#SBATCH--job-name=extracReadFromFasta_MaizeWang2020.%A_%a
#SBATCH--output=/blue/mcintyre/share/transcript_distance/analysis/SLURM_LOGS/extracReadFromFasta_MaizeWang2020.%A_%a.out
#SBATCH--error=/blue/mcintyre/share/transcript_distance/analysis/SLURM_LOGS/extracReadFromFasta_MaizeWang2020.%A_%a.err
#SBATCH --time=01:00:00
#SBATCH --ntasks=1
#SBATCH --mem=4gb
#SBATCH --array=1-12

module load samtools/1.10
module load python/3.6.5

PROJ=/blue/mcintyre/share/transcript_distance/analysis/MaizeWang2020
INDLIST=$PROJ/baseline_read_filter
INDCLS=$PROJ/SQANTI3_QC_output
PATHSCRIPT=/blue/mcintyre/share/transcript_distance/scripts
OUTD=$PROJ/SQANTI3_filter
    mkdir -p $OUTD

translist=$(ls $INDLIST/*_read_pass_filter.txt | awk -v row="$SLURM_ARRAY_TASK_ID" 'NR==row {print}')
name=$(basename $translist "_read_pass_filter.txt")

python ${PATHSCRIPT}/Common_extracFeatureFromList.py --list $translist --type fasta --input $INDCLS/${name}/${name}_corrected.fasta --out $OUTD/$name.passed.fasta
