#!/bin/bash
#SBATCH --mail-type=FAIL,END
#SBATCH--job-name=SQANTI_filter_MaizeWang.%A_%a
#SBATCH--output=/blue/mcintyre/share/transcript_distance/analysis/SLURM_LOGS/SQANTI_filter_MaizeWang.%A_%a.out
#SBATCH--error=/blue/mcintyre/share/transcript_distance/analysis/SLURM_LOGS/SQANTI_filter_MaizeWang.%A_%a.err 
#SBATCH --time=00:30:00
#SBATCH --ntasks=1
#SBATCH --mem=8gb



PROJ=/blue/mcintyre/share/transcript_distance/analysis/MaizeWang2020
cls=$PROJ/SQANTI3_filter/B73_EM.merged.passed_classification.txt
fa=$PROJ/SQANTI3_filter/B73_EM.merged.passed.fasta
gtf=$PROJ/baseline_read_filter/B73_EM.merged_corrected_associated_gene_read_filter.gtf
TOOLDIR=/blue/mcintyre/share/transcript_distance/tool/SQANTI3
CONDACONFIG=~/.conda/envs

module load conda/4.8.4
source $HPC_CONDA_BIN/activate $CONDACONFIG/SQANTI3.env
python $TOOLDIR/sqanti3_RulesFilter.py $cls $fa $gtf
