#!/bin/bash
#SBATCH --mail-user=zihaoliu@ufl.edu
#SBATCH --mail-type=FAIL,END
#SBATCH--job-name=summary_limaResult.%A_%a
#SBATCH--output=/blue/mcintyre/share/transcript_distance/analysis/SLURM_LOGS/summary_limaResult.%A_%a.out
#SBATCH--error=/blue/mcintyre/share/transcript_distance/analysis/SLURM_LOGS/summary_limaResult.%A_%a.err 
#SBATCH --cpus-per-task=12
#SBATCH --time=99:00:00
#SBATCH --ntasks=1
#SBATCH --mem=6gb


PROJ=/blue/mcintyre/share/transcript_distance/analysis/MaizeWang2020/lima
PATHSCRIPT=/blue/mcintyre/share/transcript_distance/scripts
OUTD=$PROJ/summary
    mkdir -p $OUTD

module load python/3.6.5
python ${PATHSCRIPT}/summary_lima_MaizeWang2020.py --dir $PROJ --out $OUTD --preffix MaizeWang2020
