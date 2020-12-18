#!/bin/bash
#SBATCH --mail-user=zihaoliu@ufl.edu
#SBATCH --mail-type=FAIL,END
#SBATCH--job-name=runSQANTIfilter_MaizeWang2020.%A_%a
#SBATCH--output=/blue/mcintyre/share/transcript_distance/analysis/SLURM_LOGS/runSQANTIfilter_MaizeWang2020.%A_%a.out
#SBATCH--error=/blue/mcintyre/share/transcript_distance/analysis/SLURM_LOGS/runSQANTIfilter_MaizeWang2020.%A_%a.err
#SBATCH --time=01:00:00
#SBATCH --ntasks=1
#SBATCH --mem=16gb
#SBATCH --array=1-12





PROJ=/blue/mcintyre/share/transcript_distance/analysis/MaizeWang2020
TOOLDIR=/blue/mcintyre/share/transcript_distance/tool/SQANTI3
CONDACONFIG=/home/zihaoliu/.conda/envs
cls=$(ls $PROJ/SQANTI3_filter/*.passed_classification.txt | awk -v row="$SLURM_ARRAY_TASK_ID" 'NR==row {print}')
name=$(basename $cls ".passed_classification.txt")
fa=$PROJ/SQANTI3_filter/${name}.passed.fasta
gtf=$PROJ/baseline_read_filter/${name}_corrected_associated_gene_read_filter.gtf


module load conda/4.8.4
TOOLDIR=/blue/mcintyre/share/transcript_distance/tool/SQANTI3
CONDACONFIG=/home/zihaoliu/.conda/envs

source $HPC_CONDA_BIN/activate $CONDACONFIG/SQANTI3.env
python $TOOLDIR/sqanti3_RulesFilter.py $cls $fa $gtf
