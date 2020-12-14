#!/bin/bash
#SBATCH --mail-type=FAIL,END
#SBATCH--job-name=gff2gtfAnno.%A_%a
#SBATCH--output=/blue/mcintyre/share/transcript_distance/analysis/SLURM_LOGS/gff2gtfAnno.%A_%a.out
#SBATCH--error=/blue/mcintyre/share/transcript_distance/analysis/SLURM_LOGS/gff2gtfAnno.%A_%a.err 
#SBATCH --time=1:00:00
#SBATCH --ntasks=1
#SBATCH --mem=3gb


module load gffread/0.11.8




PROJCelegans=/blue/mcintyre/share/references/Celegans_WBcel235/FSM_consolidation_Celegans_WBcel235
gffCelegans=$PROJCelegans/Celegans_WBcel235_FSM_consolidation_all_transcript.gff
gtfCelegans=$(basename $gffCelegans ".gff").gtf
sort -k1,1 -k3,3 -k4,4n -k5,5nr ${gffCelegans} | gffread - -T -o $PROJCelegans/$gtfCelegans

#PROJMaize=/blue/mcintyre/share/references/maize_b73/FSM_consolidation_maize_B73
#gffMaize=$PROJMaize/maize_B73_FSM_consolidation_all_transcript.gff
#gtfMaize=$(basename $gffMaize ".gff").gtf
#sort -k1,1 -k3,3 -k4,4n -k5,5nr ${gffMaize} | gffread - -T -o ${gtfMaize}






