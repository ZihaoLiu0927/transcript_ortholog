#!/bin/sh
#SBATCH --mail-user=adalena.nanni@ufl.edu
#SBATCH --account=mcintyre
#SBATCH --qos=mcintyre
#SBATCH --job-name=iso_format_gff
#SBATCH --mail-type=FAIL,END
#SBATCH --output=/blue/mcintyre/share/transcript_distance/analysis/SLURM_LOGS/create_formatted_gff3_isoseq3_cluster_SQANTI3_QC_GTF_combined_%A_%a.out
#SBATCH --time=12:00:00
#SBATCH --ntasks=1              # single cpu core
#SBATCH --cpus-per-task=1
#SBATCH --mem=4gb               # per processor memory
#SBATCH --array=1-2

module load conda

## Activate conda environment
## Environment created using conda env create --file /blue/mcintyre/share/transcript_distance/scripts/event_analysis/event_analysis_updated.yaml
source activate event_analysis_updated

### isoseq3 cluster SQANTI3 QC output C. elegans and Maize B73 shell script to convert GTF/GFF to Flybase-style GFF3

## Set project directory
    PROJ=/blue/mcintyre/share/transcript_distance

## Design file

    DESIGN_FILE=$PROJ/design_files/design_Celegans_MaizeWang2020_project_only_noHeader.csv

    DESIGN=$(cat $DESIGN_FILE | head -n ${SLURM_ARRAY_TASK_ID} | tail -n 1)
    IFS=',' read -ra ARRAY <<< "$DESIGN"

    DATA=${ARRAY[0]}

## Set input and output directories
    INDIR=$PROJ/analysis/${DATA}/SQANTI3_QC_isoseq3Cluster/${DATA}
    OUTDIR=$PROJ/analysis/${DATA}/EA_output_isoseq3Cluster_combined
    mkdir -p ${OUTDIR}
    echo -e "GTF/GFF formatting for ${DATA}...\nOutput in ${OUTDIR}"

### Path to GTF/GFF file to convert
    GFF=${INDIR}/${DATA}_isoseq3Cluster_corrected_associated_gene.gtf

### Output GFF path and name. Make sure this ends with ".gff"
    GFFOUT=${OUTDIR}/${DATA}_isoseq3Cluster_corrected_associated_gene.converted.gff

### Path to Event Analysis install. This is the folder containing "docs", "src", run_buildAnnotations.sh, etc.
    EVENTDIR=$PROJ/scripts/event_analysis
    SCRIPTS=${EVENTDIR}/src

### Set temp directory
    TMPDIR=${OUTDIR}/roz_${DATA}_gtf
    if [[ ! -e ${TMPDIR} ]]; then mkdir -p ${TMPDIR}; fi

#### Checking if a gff.db file exists, and if not then create one
    echo "Checking if user-supplied GTF/GFF file has a pre-generated database file"
        if [ ! -e ${GFF}.db ]
        then
        echo "No database file found! Generating..."
        python $SCRIPTS/make_gff_db_02avn.py --gff $GFF
        echo "Database generated!"
        else
        echo "Database file found! Skipping generation."
        fi

### Convert user-supplied GTF/GFF file into GFF3 format
    echo "Converting user-supplied GTF/GFF to GFF3 format"

    if [ ! -e ${GFFOUT} ]
    then
    python $SCRIPTS/convertGTF2GFF3.py --input $GFF --output ${GFFOUT}
    sort -k1,1 -k4n -k5n ${GFFOUT} > ${TMPDIR}/temp.gff
    mv ${TMPDIR}/temp.gff ${GFFOUT}
    else
    echo "Converted file found! Skipping conversion."
    fi

    if [ ! -e ${GFFOUT}.db ]
    then
    python $SCRIPTS/make_gff_db_02avn.py --gff ${GFFOUT}
    else
    echo "Converted database file found! Skipping generation."
    fi

    rm -r ${TMPDIR}

