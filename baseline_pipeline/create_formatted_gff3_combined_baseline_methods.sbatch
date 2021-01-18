#!/bin/sh
#SBATCH --mail-user=adalena.nanni@ufl.edu
#SBATCH --account=mcintyre
#SBATCH --qos=mcintyre-b
#SBATCH --job-name=comb_format_gff
#SBATCH --mail-type=FAIL,END
#SBATCH --output=/blue/mcintyre/share/transcript_distance/analysis/SLURM_LOGS/create_formatted_gff3_combined_%A_%a.out
#SBATCH --time=24:00:00
#SBATCH --ntasks=1              # single cpu core
#SBATCH --cpus-per-task=1
#SBATCH --mem=8gb               # per processor memory
#SBATCH --array=3,5,6,8-19

module load conda

## Activate conda environment
## Environment created using conda env create --file /blue/mcintyre/share/transcript_distance/scripts/event_analysis/event_analysis_updated.yaml
source activate event_analysis_updated

### Combined baseline and method output C. elegans and Maize B73 shell script to convert GTF/GFF to Flybase-style GFF3

## Set project directory
    PROJ=/blue/mcintyre/share/transcript_distance

## Design file

    DESIGN_FILE=$PROJ/design_files/design_Celegans_MaizeWang2020_noHeader.csv

    DESIGN=$(cat $DESIGN_FILE | head -n ${SLURM_ARRAY_TASK_ID} | tail -n 1)
    IFS=',' read -ra ARRAY <<< "$DESIGN"

    DATA=${ARRAY[0]}
    SAMPLE=${ARRAY[1]}
    NAME=${ARRAY[2]}

    ## Fix NAME values for MaizeWang2020
    ## FLAIR output does not contain the ".merged"
    if [[ ${DATA} == "MaizeWang2020" ]]; then
        NAME=$(basename ${NAME} .merged)
    fi

### Set input and output directories
    INDIR=$PROJ/analysis/${DATA}
    OUTDIR=$PROJ/analysis/${DATA}/EA_output_combine_baseline_methods
    ROZ=${OUTDIR}/roz_gff_concat
    mkdir -p ${ROZ}
    echo -e "GTF/GFF formatting for ${DATA}...${NAME}\nOutput in ${OUTDIR}"

### Add prefix to baseline and method GFF transcript_id values
    for TYPE in baseline FLAIR isoseq3Cluster SQANTI3filter; do
        if [[ ${TYPE} == 'baseline' ]]; then
            FILE=${INDIR}/EA_output_baseline/${NAME}*_corrected_associated_gene_read_filter.converted.gff
        elif [[ ${TYPE} == 'FLAIR' ]]; then
            FILE=${INDIR}/EA_output_FLAIR/${NAME}_flair.isoforms.converted.gff
        elif [[ ${TYPE} == 'isoseq3Cluster' ]]; then
            FILE=${INDIR}/EA_output_isoseq3Cluster/${NAME}_isoseq3Cluster_corrected_associated_gene.converted.gff
        else
            FILE=${INDIR}/EA_output_SQANTI3_filter/${NAME}*.passed_classification.filtered_lite.converted.gff
        fi
        ## Add prefix to exon and transcript elements
        ## (gene and chromosome elements dropped)
        awk -v prefix=${TYPE} 'BEGIN{OFS="\t"}$1!="##gff-version"{ \
            if($3=="exon") \
            { \
                split($9,attr,";Parent="); \
                split(attr[2],attrLess,";"); \
                num=split(attrLess[1],xcrpts,","); \
                for(count=1; count<=num; ++count) \
                { \
                    if(count==1) \
                    { \
                        newAttr=attr[1]";Parent="prefix"_"xcrpts[count]; \
                    } \
                    else \
                    { \
                        newAttr=newAttr","prefix"_"xcrpts[count]; \
                    } \
                } \
                $9=newAttr";"attrLess[2]; \
                print $0; \
            } \
            else if($3=="transcript") \
            { \
                split($9,attr,"="); \
                $9=attr[1]"="prefix"_"attr[2]"="prefix"_"attr[3]"="attr[4]; \
                print $0; \
            } \
        }' ${FILE} > ${ROZ}/${NAME}_${TYPE}.gff
    done

### Concatenate baseline, FLAIR, isoseq3 cluster, and SQANTI3 fliter GFF files
    echo "Concatenate baseline, FLAIR, isoseq3 cluster, and SQANTI3 fliter GFF files"
    cat ${ROZ}/${NAME}_*.gff | sort -k1,1 -k4,4n | uniq  \
        > ${OUTDIR}/${NAME}_combine_baseline_FLAIR_isoseq3Cluster_SQANTI3filter.gff

### Path to GTF/GFF file to convert
    GFF=${OUTDIR}/${NAME}_combine_baseline_FLAIR_isoseq3Cluster_SQANTI3filter.gff

### Output GFF path and name. Make sure this ends with ".gff"
    GFFOUT=${OUTDIR}/${NAME}_combine_baseline_FLAIR_isoseq3Cluster_SQANTI3filter.converted.gff

### Path to Event Analysis install. This is the folder containing "docs", "src", run_buildAnnotations.sh, etc.
    EVENTDIR=$PROJ/scripts/event_analysis
    SCRIPTS=${EVENTDIR}/src

### Set temp directory
    TMPDIR=${OUTDIR}/roz_${NAME}_gtf
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

