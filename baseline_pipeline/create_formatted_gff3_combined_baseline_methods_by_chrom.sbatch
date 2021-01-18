#!/bin/sh
#SBATCH --mail-user=adalena.nanni@ufl.edu
#SBATCH --account=mcintyre
#SBATCH --qos=mcintyre
#SBATCH --job-name=comb_format_gff
#SBATCH --mail-type=FAIL,END
#SBATCH --output=/blue/mcintyre/share/transcript_distance/analysis/SLURM_LOGS/create_formatted_gff3_combined_%A_%a.out
#SBATCH --time=12:00:00
#SBATCH --ntasks=1              # single cpu core
#SBATCH --cpus-per-task=1
#SBATCH --mem=4gb               # per processor memory
#SBATCH --array=1-15,23-1090

## first test: 16-22 is for Celegans,L2.gtf,L2

module load conda

## Activate conda environment
## Environment created using conda env create --file /blue/mcintyre/share/transcript_distance/scripts/event_analysis/event_analysis_updated.yaml
source activate event_analysis_updated

### Combined baseline and method output C. elegans and Maize B73 shell script to convert GTF/GFF to Flybase-style GFF3

## Set project directory
    PROJ=/blue/mcintyre/share/transcript_distance

## Design file

    DESIGN_FILE=$PROJ/design_files/design_Celegans_MaizeWang2020_chrom_noHeader.csv

    DESIGN=$(cat $DESIGN_FILE | head -n ${SLURM_ARRAY_TASK_ID} | tail -n 1)
    IFS=',' read -ra ARRAY <<< "$DESIGN"

    DATA=${ARRAY[0]}
    SAMPLE=${ARRAY[1]}
    NAME=${ARRAY[2]}
    CHROM=${ARRAY[3]}

    ## Fix NAME values for MaizeWang2020
    ## FLAIR output does not contain the ".merged"
    if [[ ${DATA} == "MaizeWang2020" ]]; then
        NAME=$(basename ${NAME} .merged)
    fi

### Set input and output directories
    INDIR=$PROJ/analysis/${DATA}
    OUTDIR=$PROJ/analysis/${DATA}/EA_output_combine_baseline_methods_by_chrom
    ROZ=${OUTDIR}/roz_gff_concat
    mkdir -p ${ROZ}
    echo -e "GTF/GFF formatting for ${DATA}...${NAME} ${CHROM}\nOutput in ${OUTDIR}"

### Add prefix to baseline and method GFF transcript_id values
    for TYPE in baseline FLAIR isoseq3Cluster SQANTI3filter; do
        if [[ ${TYPE} == 'baseline' ]]; then
            FILE=${INDIR}/sam2gtf/${NAME}.gtf
        elif [[ ${TYPE} == 'FLAIR' ]]; then
            FILE=${INDIR}/FLAIR/${NAME}_flair.isoforms.gtf
        elif [[ ${TYPE} == 'isoseq3Cluster' ]]; then
            FILE=${INDIR}/SQANTI3_QC_isoseq3Cluster/${NAME}/${NAME}_isoseq3Cluster_corrected_associated_gene.gtf
        else
            FILE=${INDIR}/SQANTI3_filter/${NAME}*.passed_classification.filtered_lite.gtf
        fi
        ## Add prefix to exon and transcript elements
        ## (gene and chromosome elements dropped)
        awk -v prefix=${TYPE} -v chrom=${CHROM} '$1==chrom \
            { \
                if($3=="exon" || $3=="transcript") \
                { \
                    split($12,xcrpt,"\""); \
                    $12="\""prefix"_"xcrpt[2]"\";"; \
                    print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9" "$10" "$11" "$12 \
                } \
            }' ${FILE} > ${ROZ}/${NAME}_${CHROM}_${TYPE}.gtf
    done

### Concatenate baseline, FLAIR, isoseq3 cluster, and SQANTI3 fliter GFF files
    echo "Concatenate baseline, FLAIR, isoseq3 cluster, and SQANTI3 fliter GFF files"
    cat ${ROZ}/${NAME}_${CHROM}_*.gtf \
        > ${OUTDIR}/${NAME}_${CHROM}_combine_baseline_FLAIR_isoseq3Cluster_SQANTI3filter.gtf

### Path to GTF/GFF file to convert
    GFF=${OUTDIR}/${NAME}_${CHROM}_combine_baseline_FLAIR_isoseq3Cluster_SQANTI3filter.gtf

### Output GFF path and name. Make sure this ends with ".gff"
    GFFOUT=${OUTDIR}/${NAME}_${CHROM}_combine_baseline_FLAIR_isoseq3Cluster_SQANTI3filter.converted.gff

### Path to Event Analysis install. This is the folder containing "docs", "src", run_buildAnnotations.sh, etc.
    EVENTDIR=$PROJ/scripts/event_analysis
    SCRIPTS=${EVENTDIR}/src

### Set temp directory
    TMPDIR=${OUTDIR}/roz_${NAME}_${CHROM}_gtf
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

