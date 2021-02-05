#!/bin/bash
#SBATCH --mail-user=zihaoliu@ufl.edu
#SBATCH --mail-type=FAIL,END
#SBATCH--job-name=runSQANTIfilter.%A_%a
#SBATCH--output=/blue/mcintyre/share/transcript_distance/analysis/SLURM_LOGS/runSQANTIfilter.%A_%a.out
#SBATCH--error=/blue/mcintyre/share/transcript_distance/analysis/SLURM_LOGS/runSQANTIfilter.%A_%a.err
#SBATCH --time=10:00:00
#SBATCH --ntasks=1
#SBATCH --mem=16gb
#SBATCH --array=1-21



PROJ=/blue/mcintyre/share/transcript_distance
DESIGN_FILE=$PROJ/design_files/design_Celegans_MaizeWang2020_noHeader.csv
DESIGN=$(cat $DESIGN_FILE | head -n $SLURM_TASK_ID | tail -n 1)
IFS=',' read -ra ARRAY <<< "$DESIGN"
SPECIES=${ARRAY[0]}
NAME=${ARRAY[2]}

TOOLDIR=$PROJ/scripts/SQANTI3
CONDACONFIG=/home/zihaoliu/.conda/envs

cls=$PROJ/analysis/$SPECIES/SQANTI3_filter/${NAME}.passed_classification.txt
fa=$PROJ/analysis/$SPECIES/SQANTI3_filter/${NAME}.passed.fasta
gtf=$PROJ/analysis/$SPECIES/baseline_read_filter/${NAME}_corrected_associated_gene_read_filter.gtf

module load conda/4.8.4

source $HPC_CONDA_BIN/activate $CONDACONFIG/SQANTI3.env
python $TOOLDIR/sqanti3_RulesFilter.py $cls $fa $gtf




