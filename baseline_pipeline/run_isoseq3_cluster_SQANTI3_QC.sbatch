#!/bin/sh
#SBATCH --mail-user=adalena.nanni@ufl.edu
#SBATCH --account=mcintyre
#SBATCH --qos=mcintyre
#SBATCH --job-name=iso_SQANTI_QC
#SBATCH --mail-type=FAIL,END
#SBATCH --output=/blue/mcintyre/share/transcript_distance/analysis/SLURM_LOGS/isoseq3_cluster_SQANTI_QC_%A_%a.out
#SBATCH --time=36:00:00
#SBATCH --ntasks=1              # single cpu core
#SBATCH --cpus-per-task=1
#SBATCH --mem=4gb               # per processor memory
#SBATCH --array=1-19

module load conda

### Activate conda environment
### Environment created using conda env create --file /blue/mcintyre/share/transcript_distance/scripts/SQANTI3.conda_env.yaml
source activate SQANTI3.env

### SQANTI3 QC of isoseq3 cluster output C. elegans and Maize B73 shell script

### Set project directory
    PROJ=/blue/mcintyre/share/transcript_distance

### Set up python path
    export PYTHONPATH=$PYTHONPATH:$PROJ/scripts/cDNA_Cupcake/sequence/
    export PYTHONPATH=$PYTHONPATH:$PROJ/scripts/cDNA_Cupcake/

### Design file

    DESIGN_FILE=$PROJ/design_files/design_Celegans_MaizeWang2020_noHeader.csv

    DESIGN=$(cat $DESIGN_FILE | head -n ${SLURM_ARRAY_TASK_ID} | tail -n 1)
    IFS=',' read -ra ARRAY <<< "$DESIGN"

    DATA=${ARRAY[0]}
    SAMPLE=${ARRAY[1]}
    NAME=${ARRAY[2]}

### Set input and output directories
    INDIR=$PROJ/analysis/${DATA}/isoseq3cluster/isoseqClusterGTF
    OUTDIR=$PROJ/analysis/${DATA}/SQANTI3_QC_isoseq3cluster
    mkdir -p ${OUTDIR}
    echo -e "SQANTI3 QC of ${DATA}...${NAME}\nOutput in ${OUTDIR}"

### Set proper FSM consolidated reference GTF and genome FASTA
    if [[ ${DATA} == Celegans ]]; then
        REF=/blue/mcintyre/share/references/Celegans_WBcel235/FSM_consolidation_Celegans_WBcel235/Celegans_WBcel235_FSM_consolidation_all_transcript.gtf
        FA=/blue/mcintyre/share/references/Celegans_WBcel235/Caenorhabditis_elegans.WBcel235.dna.toplevel.fa
    else
        REF=/blue/mcintyre/share/references/maize_b73/FSM_consolidation_maize_B73/maize_B73_FSM_consolidation_all_transcript.gtf
        FA=/blue/mcintyre/share/references/maize_b73/b73_genome_ensemblgenomes/fasta/zea_mays/dna/Zea_mays.B73_RefGen_v4.dna.toplevel.fa
    fi

### Path to input isoseq3 cluster GTF
    ISO=${INDIR}/${NAME}.isoseq3cluster_corrected.gtf

### Check if SQANTI3 report has already been generated
    if [[ -e ${OUTDIR}/${NAME}/${NAME}_isoseq3Cluster_sqanti_report.pdf ]]; then
        echo "!!! ${OUTDIR}/${NAME}/${NAME}_isoseq3Cluster_sqanti_report.pdf already generated"
        continue
    else
        ## SQANTI3 QC
        ##     --gtf option indicates GTF isoform input
        ##     --skipORF option to skip ORF prediction (runs much faster)
        ##     --force_id_ignore option to allow transcript ID names that are not like PB.X.Y
        ## Order of files is isoforms (FA/FQ/GTF), annotation (ref GTF/GFF), genome(ref FA)
        python $PROJ/scripts/SQANTI3/sqanti3_qc.py \
            --gtf \
            --skipORF \
            --force_id_ignore \
            --dir ${OUTDIR}/${NAME} \
            -o ${NAME}_isoseq3Cluster \
            ${ISO} \
            ${REF} \
            ${FA}
    fi

### Check that SQANTI3 GTF and classification files were generated
    if [[ ! -e ${OUTDIR}/${NAME}/${NAME}_isoseq3Cluster_classification.txt || ! -e ${OUTDIR}/${NAME}/${NAME}_isoseq3Cluster_corrected.gtf ]]; then
        echo "!!! ERROR: REQUIRED SQANTI3 QC OUTPUT FILES DO NOT EXIST"
        continue
### Check if gene_id corrected GTF has already been generated
    elif [[ -e ${OUTDIR}/${NAME}/${NAME}_isoseq3Cluster_corrected_associated_gene.gtf ]]; then
        echo "!!! ${NAME}_isoseq3Cluster_corrected_associated_gene.gtf already generated"
        continue
    else
        ## Get associated gene values for each isoform from classification files
        ## Correct the gene_id values of the GTF output from SQANTI3 QC to
        ##     match the gene_id values associated with the transcripts in
        ##     the classification file
        python $PROJ/scripts/correct_SQANTI_QC_GTF_gene_id.py \
            -c ${OUTDIR}/${NAME}/${NAME}_isoseq3Cluster_classification.txt \
            -g ${OUTDIR}/${NAME}/${NAME}_isoseq3Cluster_corrected.gtf \
            --output-classification ${OUTDIR}/${NAME}_isoseq3Cluster_no_associated_gene_classification.txt \
            --output-gtf ${OUTDIR}/${NAME}_isoseq3Cluster_no_associated_gene.gtf \
            -o ${OUTDIR}/${NAME}/${NAME}_isoseq3Cluster_corrected_associated_gene.gtf
        echo -e "\t...${NAME} gene_id corrected GTF output to ${NAME}_isoseq3Cluster_corrected_associated_gene.gtf"
        if [[ $(wc -l ${OUTDIR}/${NAME}_isoseq3Cluster_no_associated_gene.gtf | awk '{print $1}') == 0 ]]; then
            rm ${OUTDIR}/${NAME}_isoseq3Cluster_no_associated_gene_classification.txt
            rm ${OUTDIR}/${NAME}_isoseq3Cluster_no_associated_gene.gtf
        fi
    fi

