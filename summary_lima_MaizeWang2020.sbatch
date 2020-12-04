#!/bin/bash

PROJ=/blue/mcintyre/share/transcript_distance/analysis/MaizeWang2020/lima
PATHSCRIPT=/blue/mcintyre/share/transcript_distance/scripts

module load python/3.6.5
python ${PATHSCRIPT}/summary_lima_MaizeWang2020.py --dir $PROJ --out $PROJ/summary --preffix MaizeWang2020
