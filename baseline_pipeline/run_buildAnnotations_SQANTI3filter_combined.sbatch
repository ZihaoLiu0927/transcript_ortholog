#!/bin/sh
#SBATCH --mail-user=adalena.nanni@ufl.edu
#SBATCH --account=mcintyre
#SBATCH --qos=mcintyre-b
#SBATCH --job-name=SQANTI_build_annot
#SBATCH --mail-type=FAIL,END
#SBATCH --output=/blue/mcintyre/share/transcript_distance/analysis/SLURM_LOGS/buildAnnotations_SQANTI3filter_combined_%A_%a.out
#SBATCH --time=26:00:00
#SBATCH --ntasks=1              # single cpu core
#SBATCH --cpus-per-task=1
#SBATCH --mem=48gb               # per processor memory
#SBATCH --array=1-2

module load conda

## Activate conda environment
## Environment created using conda env create --file /blue/mcintyre/share/transcript_distance/scripts/event_analysis/event_analysis_updated.yaml
source activate event_analysis_updated

### SQANTI3 filter output C. elegans and Maize B73 shell script to build Event Analysis annotations

## Set project directory
    PROJ=/blue/mcintyre/share/transcript_distance

## Design file

    DESIGN_FILE=$PROJ/design_files/design_Celegans_MaizeWang2020_project_only_noHeader.csv

    DESIGN=$(cat $DESIGN_FILE | head -n ${SLURM_ARRAY_TASK_ID} | tail -n 1)
    IFS=',' read -ra ARRAY <<< "$DESIGN"

    DATA=${ARRAY[0]}

    echo -e "EA annotation building for ${DATA}...\nOutput in ${OUTDIR}"

### Path to Event Analysis install
    EVENTDIR=$PROJ/scripts/event_analysis

### Label to prefix to annotation files (e.g. mm10_refseq, dmel617, hg38_ens)
    PREFIX=${DATA}_SQANTI3filter

### Output directory. If it does not exist, the annotation build script with create it
    OUTDIR=$PROJ/analysis/${DATA}/EA_output_SQANTI3_filter_combined/${DATA}_SQANTI3filter_150bp_annotations

### Path to formatted GFF3 file (should in a FlyBase GFF3 format)
    GFF=$PROJ/analysis/${DATA}/EA_output_SQANTI3_filter_combined/${DATA}.passed_classification.filtered_lite.converted.gff

### Path to genome FASTA
    if [[ ${DATA} == "Celegans" ]]; then
        FASTA=/blue/mcintyre/share/references/Celegans_WBcel235/Caenorhabditis_elegans.WBcel235.dna.toplevel.fa
    elif [[ ${DATA} == "MaizeWang2020" ]]; then
        FASTA=/blue/mcintyre/share/references/maize_b73/b73_genome_ensemblgenomes/fasta/zea_mays/dna/Zea_mays.B73_RefGen_v4.dna.toplevel.fa
    fi

### Size (in nt) of reads or each read pair (set to maximum if you have reads of
### various lengths, i.e. from adapter trimming)
    READSIZE=150

### Build Event Analysis annotations
    cd ${EVENTDIR}
    sh ./run_buildAnnotations.sh ${PREFIX} ${GFF} ${FASTA} ${READSIZE} ${OUTDIR}

### Remove the fasta files of EA output to reduce memory
rm ${OUTDIR}/*.fa
