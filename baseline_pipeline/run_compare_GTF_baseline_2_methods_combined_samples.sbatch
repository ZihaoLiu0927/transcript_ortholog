#!/bin/sh
#SBATCH --mail-user=adalena.nanni@ufl.edu
#SBATCH --account=mcintyre
#SBATCH --qos=mcintyre
#SBATCH --job-name=compare_base2methods
#SBATCH --mail-type=FAIL,END
#SBATCH --output=/blue/mcintyre/share/transcript_distance/origin_analysis/SLURM_LOGS/compare_baseline_2_methods_combined_%A_%a.out
#SBATCH --time=03:00:00
#SBATCH --ntasks=1              # single cpu core
#SBATCH --cpus-per-task=1
#SBATCH --mem=150gb               # per processor memory
#SBATCH --array=1
##SBATCH --partition=bigmem

module load python/3.7

## Set read threshold
## Number of reads in gene for evidence
NUM=10

## Get lists of reads present in baseline and each method
##     (SQANTI3 filter, isoseq3 cluster, and FLAIR)
##     FLAIR and isoseq3 replace read ID value with reference transcript_id
##     values (FLAIR) or new clustered transcript_id values (isoseq3 cluster)
## Count how many reads are retained or modified in each method and across methods

## Set project directory
PROJ=/blue/mcintyre/share/transcript_distance

## Set summary file
SUMMARY=$PROJ/origin_analysis/compare_baseline_to_methods_summary_counts_combined_samples.csv

echo "dataset,sample,sum_in_baseline_read,sum_in_SQANTI3filter,sum_in_isoseq3Cluster,sum_in_FLAIR,sum_in_SQANTI3filter_isoseq3Cluster,sum_in_SQANTI3filter_FLAIR,sum_in_isoseq3Cluster_FLAIR,sum_in_SQANTI3filter_isoseq3Cluster_FLAIR,sum_in_SQANTI3filter_not_in_FLAIR,sum_in_SQANTI3filter_not_in_isoseq3Cluster,sum_in_FLAIR_not_in_SQANTI3filter,sum_in_FLAIR_not_in_isoseq3Cluster,sum_in_isoseq3Cluster_not_in_SQANTI3filter,sum_in_isoseq3Cluster_not_in_FLAIR,sum_modified_in_isoseq3Cluster,sum_modified_internal_in_isoseq3Cluster,sum_modified_in_FLAIR,sum_modified_internal_in_FLAIR,sum_modified_in_isoseq3Cluster_FLAIR,sum_modified_internal_in_isoseq3Cluster_FLAIR,num_baseline_gene,num_single_read_baseline_gene,num_gene_at_least_${NUM}_reads_per_gene,num_gene_less_than_${NUM}_reads_per_gene,num_at_least_${NUM}_read_baseline_gene_in_SQANTI3filter,num_at_least_${NUM}_read_baseline_gene_in_isoseq3Cluster_no_mutligene_cluster,num_at_least_${NUM}_read_baseline_gene_in_FLAIR_no_mutligene_cluster,perc_at_least_${NUM}_read_baseline_gene_in_SQANTI3filter,perc_at_least_${NUM}_read_baseline_gene_in_isoseq3Cluster_no_mutligene_cluster,perc_at_least_${NUM}_read_baseline_gene_in_FLAIR_no_mutligene_cluster,num_multigene_cluster_in_isoseq3Cluster,num_multigene_cluster_in_FLAIR,num_single_gene_cluster_in_isoseq3Cluster,num_single_gene_cluster_in_FLAIR,num_multi_read_baseline_gene,num_multi_read_baseline_gene_in_SQANTI3filter,num_multi_read_baseline_gene_in_isoseq3Cluster_no_mutligene_cluster,num_multi_read_baseline_gene_in_FLAIR_no_mutligene_cluster,perc_multi_read_baseline_gene_in_SQANTI3filter,perc_multi_read_baseline_gene_in_isoseq3Cluster_no_mutligene_cluster,perc_multi_read_baseline_gene_in_FLAIR_no_mutligene_cluster" > ${SUMMARY}

## Loop over datasets
for DATA in MaizeWang2020 Celegans; do
    echo "${DATA} comparisons..."
    ## Set output directory
    OUTDIR=$PROJ/origin_analysis/${DATA}/compare_baseline_to_methods_combined
    mkdir -p ${OUTDIR}

    ## Get combined sample read filtered GTF
    SAMPLE=$PROJ/origin_analysis/${DATA}/baseline_read_filter/${DATA}_corrected_associated_gene_read_filter.gtf

    ## SQANTI3 filter name is different for MaizeWang2020
    if [[ ${DATA} == "MaizeWang2020" ]]; then
        NAME=MaizeWang2020.merged
    else
        NAME=${DATA}
    fi

    ## Output dataset and sample names to summary file
    echo -n "${DATA},combined," >> ${SUMMARY}
    ## Compare and count
    rm ${OUTDIR}/${DATA}_temp.db
    python $PROJ/scripts/compare_GTF_baseline_2_methods.py \
        -b ${SAMPLE} \
        -i $PROJ/origin_analysis/${DATA}/SQANTI3_QC_isoseq3Cluster/${DATA}/${DATA}_isoseq3Cluster_corrected_associated_gene.gtf \
        -r $PROJ/origin_analysis/${DATA}/isoseq3cluster/${DATA}.isoseq3_clustered.cluster_report.csv \
        -s $PROJ/origin_analysis/${DATA}/SQANTI3_filter/${NAME}.passed_classification.filtered_lite.gtf \
        -f $PROJ/origin_analysis/${DATA}/FLAIR/${DATA}_flair.isoforms.gtf \
        -m $PROJ/origin_analysis/${DATA}/FLAIR/${DATA}_flair.isoform.read.map.txt \
        -n ${NUM} \
        -g ${OUTDIR}/${DATA}_at_least_${NUM}_read_gene_list.csv \
        -o ${OUTDIR}/${DATA}_read_flags.csv \
        -db ${OUTDIR}/${DATA}_temp.db
    rm ${OUTDIR}/${DATA}_temp.db

    ## Get counts of comparisons
    python $PROJ/scripts/count_compare_baseline_to_methods_04avn.py \
        -i ${OUTDIR}/${DATA}_read_flags.csv \
        -c ${OUTDIR}/${DATA}_counts.txt \
        >> ${SUMMARY}

done
