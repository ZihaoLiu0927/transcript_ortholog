#!/bin/bash
#SBATCH --mail-type=FAIL,END
#SBATCH--job-name=MaizeWang2020_gff2gtfAnno.%A_%a
#SBATCH--output=/blue/mcintyre/share/transcript_distance/analysis/SLURM_LOGS/MaizeWang2020_gff2gtfAnno.%A_%a.out
#SBATCH--error=/blue/mcintyre/share/transcript_distance/analysis/SLURM_LOGS/MaizeWang2020_gff2gtfAnno.%A_%a.err 
#SBATCH --time=1:00:00
#SBATCH --ntasks=1
#SBATCH --mem=3gb


module load gffread/0.11.8

PROJMaize=/blue/mcintyre/share/references/maize_b73/FSM_consolidation_maize_B73
gffMaize=$PROJMaize/maize_B73_FSM_consolidation_all_transcript.gff
gtfMaize=$(basename $gffMaize ".gff").gtf
sort -k1,1 -k3,3 -k4,4n -k5,5nr ${gffMaize} | gffread - -T -o ${gtfMaize}






