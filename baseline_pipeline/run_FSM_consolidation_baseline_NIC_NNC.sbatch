#!/bin/sh
#SBATCH --mail-user=adalena.nanni@ufl.edu
#SBATCH --account=mcintyre
#SBATCH --qos=mcintyre
#SBATCH --job-name=NIC_NNC_consol
#SBATCH --mail-type=FAIL,END
#SBATCH --output=/blue/mcintyre/share/transcript_distance/origin_analysis/SLURM_LOGS/FSM_consolidation_baseline_NIC_NNC_%A_%a.out
#SBATCH --time=24:00:00
#SBATCH --ntasks=1              # single cpu core
#SBATCH --cpus-per-task=1
#SBATCH --mem=32gb               # per processor memory
#SBATCH --array=1

module load python/3.7

## Shell script for FSM consolidation for baseline NIC/NNC

## Set project directory
    PROJ=/blue/mcintyre/share/transcript_distance

## Design file

    DESIGN_FILE=$PROJ/design_files/design_Celegans_MaizeWang2020_project_only_noHeader.csv

    DESIGN=$(cat $DESIGN_FILE | head -n ${SLURM_ARRAY_TASK_ID} | tail -n 1)
    IFS=',' read -ra ARRAY <<< "$DESIGN"

    DATA=${ARRAY[0]}

#####################################################################
##################### Variables for user to set #####################
#####################################################################

    ## Path to FSM consolidation scripts
    SCRIPTS=$PROJ/scripts/FSM_consolidation

    ## Label to prefix to otuput files (e.g. mm10_refseq, dmel617, hg38_ens)
    PREFIX=${DATA}_NIC_NNC

    ## Label to prefix representative transcript_id values in output
    ## No spaces or special characters (e.g. tr, transcript, etc.)
    ## Default is "tr"
    XCRPT_PREFIX=novel

    ## Output directory. If it does not exist, the script will create it
    OUTDIR=$PROJ/$PROJ/origin_analysis/${DATA}/EA_output_baseline_combined/FSM_consolidation_${PREFIX}
    echo -e "GTF/GFF formatting for ${DATA}...\nOutput in ${OUTDIR}"

    ## Path to Event Analysis annotation output
    EVENTPATH=$PROJ/origin_analysis/${DATA}/EA_output_baseline_combined/${DATA}_NIC_NNC_150bp_annotations

    ## Prefix used in Event Analysis annotation files (e.g. mm10_refseq, dmel617, hg38_ens)
    EVENTPREFIX=${DATA}_NIC_NNC

    ## Path to GFF database file generated prior to Event Analysis annotations
    GFFDB=$PROJ/origin_analysis/${DATA}/EA_output_baseline_combined/${DATA}_read_filter_NIC_NNC.converted.gff.db

#####################################################################
############# Code below does not need changed by user ##############
#####################################################################

    ## Path to CSV file of event_id to transcript_id to gene_id from
    ##     Event Analysis annotations (*_event2transcript2gene_index.csv)
    EVENT=${EVENTPATH}/${EVENTPREFIX}_event2transcript2gene_index.csv

    ## Path to CSV file of junction annotations including transcript_id
    ##     from Event Analysis annotations (*_annotated_junctions.csv)
    JUNC=${EVENTPATH}/${EVENTPREFIX}_annotated_junctions.csv

    ## Path to CSV file of fragment annotations including transcript_id
    ##     from Event Analysis annotations (*_exon_fragment_annotations.csv)
    FRAG=${EVENTPATH}/${EVENTPREFIX}_exon_fragment_annotations.csv

    ## Set XCRPT_PREFIX to default if not set by user
    if [[ ${XCRPT_PREFIX}=="" ]]; then
        XCRPT_PREFIX=tr
    fi

    ## Collapse transcriptome by junctions
    cd ${SCRIPTS}
    bash ./run_FSM_consolidation.sh ${PREFIX} ${EVENT} ${JUNC} ${FRAG} ${GFFDB} ${OUTDIR}
