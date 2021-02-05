#!/bin/bash
#################################################################################
# Run Full Splice Match (FSM) Consolidation
# This is a bash script to consolidate transcriptome by shared junctions and
#     mono-exon genes following Event Analysis annotation building step.
#
# Usage:
# Running in home directory
# bash ./run_FSM_consolidation.sh ${PREFIX} ${EVENT} ${JUNC} ${FRAG} ${GFFDB} \
#                              ${OUTDIR}
#
# PREFIX       : Prefix to append to created output files
#
# XCRPT_PREFIX : Prefix for representative transcript_id values in output
#
# EVENT        : Path to CSV file of event_id to transcript_id to gene_id from
#                Event Analysis annotations (*_event2transcript2gene_index.csv)
#
# JUNC         : Path to CSV file of junction annotations including transcript_id
#                from Event Analysis annotations (*_annotated_junctions.csv)
#
# FRAG         : Path to CSV file of fragment annotations including transcript_id
#                from Event Analysis annotations (*_exon_fragment_annotations.csv)
#
# GFFDB        : Path to GFF database file generated prior to Event Analaysis
#                annotations (*.gff.db)
#
# OUTDIR       : Output directory. This will be created if it does not exist
#
#################################################################################

## Parse command line arguments
PREFIX=${1}
XCRPT_PREFIX=${2}
EVENT=${3}
JUNC=${4}
FRAG=${5}
GFFDB=${6}
OUTDIR=${7}

## Set paths
SCRIPTS=$(pwd)
if [ ! -e ${OUTDIR} ]; then mkdir -p ${OUTDIR}; fi

## Set log file - remove and replace if already exists
LOGFILE=${OUTDIR}/${PREFIX}_FSM_consolidation_log.txt
if [[ -f ${LOGFILE} ]]; then
	rm ${LOGFILE}
	touch ${LOGFILE}
fi


############################################
#### FSM CONSOLIDATION OF TRANSCRIPTOME ####
############################################

## Get start time of process
start=$(date +%s.%N)
## Extract all unique pairs of individual transcript_id to gene_id
## This is the file junction variables will be merged on
echo "Extracting unique transcript_id gene_id pairs"
/usr/bin/time -f "\tTime elapsed (mm:ss) = %E\n\tMax Mem (kb) = %M\n" \
python ${SCRIPTS}/extract_transcript_id_gene_id_pairs.py \
    -i ${EVENT} \
    -o ${OUTDIR}/${PREFIX}_transcript_id_2_gene_id.csv

## Get frequency of transcripts per gene prior to collapse
echo "Get frequency of transcripts per gene prior to collapse"
/usr/bin/time -f "\tTime elapsed (mm:ss) = %E\n\tMax Mem (kb) = %M\n" \
python ${SCRIPTS}/get_freq_transcript_per_gene_pre_collapse.py \
    -i ${OUTDIR}/${PREFIX}_transcript_id_2_gene_id.csv \
    -o ${OUTDIR}/${PREFIX}_transcript_per_gene_freq_pre_collapse.csv \
    >> ${LOGFILE}

## FSM Consolidation and mono-exon gene collapse
echo "FSM Consolidation and mono-exon gene collapse"
/usr/bin/time -f "\tTime elapsed (mm:ss) = %E\n\tMax Mem (kb) = %M\n" \
python ${SCRIPTS}/FSM_consolidation_03avn.py \
    -i ${OUTDIR}/${PREFIX}_transcript_id_2_gene_id.csv \
    -f ${FRAG} \
    -j ${JUNC} \
    -t ${XCRPT_PREFIX} \
    -d ${OUTDIR} \
    -p ${PREFIX} \
    >> ${LOGFILE}

## Output single-transcript and multi-transcript GTF for FSM consolidation
echo "Output single-transcript and multi-transcript GTF for FSM consolidation"
/usr/bin/time -f "\tTime elapsed (mm:ss) = %E\n\tMax Mem (kb) = %M\n" \
python ${SCRIPTS}/get_FSM_consolidation_GFF.py \
    -i ${OUTDIR}/${PREFIX}_FSM_consolidation.csv \
    -g ${GFFDB} \
    -d ${OUTDIR} \
    -p ${PREFIX} \
    >> ${LOGFILE}

## Concatenate single-transcript and multi-transcript GTF files and sort
echo "Concatenate single-transcript and multi-transcript GTF files and sort"
/usr/bin/time -f "\tTime elapsed (mm:ss) = %E\n\tMax Mem (kb) = %M\n" \
cat ${OUTDIR}/${PREFIX}_FSM_consolidation_single_transcript.gff \
    ${OUTDIR}/${PREFIX}_FSM_consolidation_multi_transcript.gff \
    | sort -k1,1 -k4,4n -k5,5n \
    > ${OUTDIR}/${PREFIX}_FSM_consolidation_all_transcript.gff

## Generate GFF database file - using script from Event Analysis code
echo "Generate GFF database file - using script from Event Analysis code"
/usr/bin/time -f "\tTime elapsed (mm:ss) = %E\n\tMax Mem (kb) = %M\n" \
python ${SCRIPTS}/make_gff_db.py \
    --gff ${OUTDIR}/${PREFIX}_FSM_consolidation_all_transcript.gff

## Get total time of script running
duration=$(echo "$(date +%s.%N) - $start" | bc)
hours=$(echo "$duration / 3600" | bc)
minutes=$(echo "($duration % 3600) / 60" | bc)
seconds=$(echo "($duration % 3600) % 60" | bc)
echo "
Total FSM Consolidation Runtime: $hours:$minutes:$seconds (hh:mm:ss)
"
