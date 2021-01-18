#!/bin/bash
#SBATCH --mail-type=FAIL,END
#SBATCH--job-name=MaizeB73_baseline_visualHistgram.%A_%a
#SBATCH--output=/blue/mcintyre/share/transcript_distance/analysis/SLURM_LOGS/MaizeB73_baseline_visualHistgram.%A_%a.out
#SBATCH--error=/blue/mcintyre/share/transcript_distance/analysis/SLURM_LOGS/MaizeB73_baseline_visualHistgram.%A_%a.err 
#SBATCH --cpus-per-task=12
#SBATCH --time=01:00:00
#SBATCH --ntasks=1
#SBATCH --mem=16gb

PROJ=/nfshome/zihaoliu/mclab/SHARE/McIntyre_Lab/transcript_distance/MaizeB73_baseline_visualization
INPUT=/nfshome/zihaoliu/mclab/SHARE/McIntyre_Lab/useful_maize_info/maize_B73_noMtPt_150bp

python $PROJ/run_visual_baseline.py \
	-ig $INPUT/maize_B73_event2transcript2gene_index.csv \
	-if $INPUT/maize_B73_exon_fragment_annotations.csv \
	-ie $INPUT/maize_B73_fusion_annotations.csv \
	-o $PROJ \
	-p maizeB73_baseline
