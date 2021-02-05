#!/bin/bash
#SBATCH --mail-type=FAIL,END
#SBATCH--job-name=shared_gff2gtfAnno.%A_%a
#SBATCH--output=/blue/mcintyre/share/transcript_distance/analysis/SLURM_LOGS/shared_gff2gtfAnno.%A_%a.out
#SBATCH--error=/blue/mcintyre/share/transcript_distance/analysis/SLURM_LOGS/shared_gff2gtfAnno.%A_%a.err 
#SBATCH --time=2:00:00
#SBATCH --ntasks=1
#SBATCH --mem=3gb

module load gffread/0.11.8


for species in Celegans_WBcel235 maize_B73
do
    PROJ=/blue/mcintyre/share/references/Celegans_WBcel235/FSM_consolidation_${species}
    gff=$PROJ/${species}_FSM_consolidation_all_transcript.gff
    gtf=$(basename $gff ".gff").gtf
    sort -k1,1 -k3,3 -k4,4n -k5,5nr ${gff} | gffread - -T -o $PROJ/$gtf
done






