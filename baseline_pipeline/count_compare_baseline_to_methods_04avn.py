#!/usr/bin/env python

import argparse
import pandas as pd

def getOptions():
    # Parse command line arguments
    parser = argparse.ArgumentParser(description="Get counts for comparisons of baseline reads to methods")

    # Input data
    parser.add_argument("-i", "--input", dest="inFile", required=True, help="CSV file of baseline reads and associated method flags and information")
    parser.add_argument("-n", "--number-reads", dest="inNum", required=False, type=int, choices=range(0,1000000), help="Number of reads to be used for flagging reads in genes with at least N reads evidence and output associated gene values for downstream analysis")

    # Output data
    parser.add_argument("--output-SQANTI", dest="outSQANTI", required=False, help="Output SQANTI3 filter only file to be used in comparisons to SQANTI classifications")
    parser.add_argument("-c", "--output-counts", dest="outCount", required=True, help="Path to output count file")

    args = parser.parse_args()
    return args

def main():
    # Get input file
    baselineIsoFlair = pd.read_csv(args.inFile, low_memory=False)

    # Output SQANTI3 filter only file
    # (to be used in comparisons to SQANTI classifications)
    if args.outSQANTI is not None:
        baselineIsoFlair[baselineIsoFlair['flag_in_SQANTI3filter']==1][[
                'readID','flag_modified_internal_in_FLAIR',
                'flag_modified_internal_in_isoseq3Cluster']].to_csv(args.outSQANTI, index=False)

### Counts
    # Open output file
    outFile = open(args.outCount,'w')
    # Number of reads/genes in baseline
    outFile.write("\n{} reads in baseline in {} annotated genes\n".format(
            len(baselineIsoFlair),baselineIsoFlair['gene_id'].nunique()))
    outFile.write("\n{} genes contain a single baseline read\n".format(baselineIsoFlair[baselineIsoFlair['flag_single_read_gene_baseline']==1]['gene_id'].nunique()))
    outFile.write("\n{0} of the baseline reads ({1:.2%}) are mono-exonic\n".format(
            baselineIsoFlair['flag_monoexon_baseline'].sum(),
            baselineIsoFlair['flag_monoexon_baseline'].sum()/len(baselineIsoFlair)))
    # Number of clusters in isoseq3 cluster and FLAIR that contain baseline reads from multiple genes
    outFile.write("\n{} multigene clusters retained by isoseq3 cluster\n".format(
            baselineIsoFlair[(baselineIsoFlair['flag_multigene_isoseq3Cluster']==1)&(baselineIsoFlair['flag_in_isoseq3Cluster']==1)]['iso_cluster_id'].nunique()))
    outFile.write("{} multigene clusters retained by FLAIR\n".format(
            baselineIsoFlair[(baselineIsoFlair['flag_multigene_FLAIR']==1)&(baselineIsoFlair['flag_in_FLAIR']==1)]['flair_transcript_id'].nunique()))
    # Number and % of genes with >1 baseline read per gene and no multi-gene clusters
    #     with at least one isoform retained in SQANTI3 filter, isoseq3 cluster, and FLAIR
    outFile.write("\n{} genes with > 1 baseline read associated\n".format(
            baselineIsoFlair[baselineIsoFlair['flag_single_read_gene_baseline']==0]['gene_id'].nunique()))
    outFile.write("{0} genes ({1:.2%}) with at least one isoform retained in SQANTI3 filter\n".format(
            baselineIsoFlair[(baselineIsoFlair['flag_in_SQANTI3filter']==1)&(baselineIsoFlair['flag_single_read_gene_baseline']==0)]['gene_id'].nunique(),
            baselineIsoFlair[(baselineIsoFlair['flag_in_SQANTI3filter']==1)&(baselineIsoFlair['flag_single_read_gene_baseline']==0)]['gene_id'].nunique()/baselineIsoFlair[baselineIsoFlair['flag_single_read_gene_baseline']==0]['gene_id'].nunique()))
    outFile.write("{0} genes ({1:.2%}) with only one isoform/cluster retained in isoseq3 cluster (excluding multi-gene clusters)\n".format(
            baselineIsoFlair[(baselineIsoFlair['flag_in_isoseq3Cluster']==1)&(baselineIsoFlair['flag_single_read_gene_baseline']==0)&(baselineIsoFlair['flag_multigene_isoseq3Cluster']==0)]['gene_id'].nunique(),
            baselineIsoFlair[(baselineIsoFlair['flag_in_isoseq3Cluster']==1)&(baselineIsoFlair['flag_single_read_gene_baseline']==0)&(baselineIsoFlair['flag_multigene_isoseq3Cluster']==0)]['gene_id'].nunique()/baselineIsoFlair[baselineIsoFlair['flag_single_read_gene_baseline']==0]['gene_id'].nunique()))
    outFile.write("{0} genes ({1:.2%}) with at least one isoform/cluster retained in FLAIR (excluding multi-gene clusters)\n".format(
            baselineIsoFlair[(baselineIsoFlair['flag_in_FLAIR']==1)&(baselineIsoFlair['flag_single_read_gene_baseline']==0)&(baselineIsoFlair['flag_multigene_FLAIR']==0)]['gene_id'].nunique(),
            baselineIsoFlair[(baselineIsoFlair['flag_in_FLAIR']==1)&(baselineIsoFlair['flag_single_read_gene_baseline']==0)&(baselineIsoFlair['flag_multigene_FLAIR']==0)]['gene_id'].nunique()/baselineIsoFlair[baselineIsoFlair['flag_single_read_gene_baseline']==0]['gene_id'].nunique()))
    outFile.write("{0} genes ({1:.2%}) with at least one isoform/cluster retained in SQANTI3 filter, isoseq3 cluster, and FLAIR (excluding multi-gene clusters)\n".format(
            baselineIsoFlair[(baselineIsoFlair['flag_in_SQANTI3filter_isoseq3Cluster_FLAIR']==1)&(baselineIsoFlair['flag_single_read_gene_baseline']==0)&(baselineIsoFlair['flag_multigene_FLAIR']==0)&(baselineIsoFlair['flag_multigene_isoseq3Cluster']==0)]['gene_id'].nunique(),
            baselineIsoFlair[(baselineIsoFlair['flag_in_SQANTI3filter_isoseq3Cluster_FLAIR']==1)&(baselineIsoFlair['flag_single_read_gene_baseline']==0)&(baselineIsoFlair['flag_multigene_FLAIR']==0)&(baselineIsoFlair['flag_multigene_isoseq3Cluster']==0)]['gene_id'].nunique()/baselineIsoFlair[baselineIsoFlair['flag_single_read_gene_baseline']==0]['gene_id'].nunique()))

    # Number and % retained in SQANTI3 filter
    outFile.write("\n{0} reads from baseline retained in SQANTI3 filter ({1:.2%} of baseline)\n".format(
            baselineIsoFlair['flag_in_SQANTI3filter'].sum(),
            baselineIsoFlair['flag_in_SQANTI3filter'].sum()/len(baselineIsoFlair)))
    # Number and % retained in isoseq3 clustering and number and % modified
    outFile.write("\n{0} reads from baseline retained in isoseq3 cluster ({1:.2%} of baseline)\n".format(
            baselineIsoFlair['flag_in_isoseq3Cluster'].sum(),
            baselineIsoFlair['flag_in_isoseq3Cluster'].sum()/len(baselineIsoFlair)))
    outFile.write("\n{0} reads from baseline modified by at least one nt (internal/5'/3') in isoseq3 cluster ({1:.2%} of baseline, {2:.2%} of reads retained in isoseq3 cluster)\n".format(
            baselineIsoFlair['flag_modified_in_isoseq3Cluster'].sum(),
            baselineIsoFlair['flag_modified_in_isoseq3Cluster'].sum()/len(baselineIsoFlair),
            baselineIsoFlair['flag_modified_in_isoseq3Cluster'].sum()/baselineIsoFlair['flag_in_isoseq3Cluster'].sum()))
    outFile.write("\n{0} reads from baseline modified internally by at least one nt in isoseq3 cluster ({1:.2%} of baseline, {2:.2%} of reads retained in isoseq3 cluster, {3:.2%} of reads modified in isoseq3 cluster)\n".format(
            baselineIsoFlair['flag_modified_internal_in_isoseq3Cluster'].sum(),
            baselineIsoFlair['flag_modified_internal_in_isoseq3Cluster'].sum()/len(baselineIsoFlair),
            baselineIsoFlair['flag_modified_internal_in_isoseq3Cluster'].sum()/baselineIsoFlair['flag_in_isoseq3Cluster'].sum(),
            baselineIsoFlair['flag_modified_internal_in_isoseq3Cluster'].sum()/baselineIsoFlair['flag_modified_in_isoseq3Cluster'].sum()))
    # Number and % retained in FLAIR and number and % modified
    outFile.write("\n{0} reads from baseline retained in FLAIR ({1:.2%} of baseline)\n".format(
            baselineIsoFlair['flag_in_FLAIR'].sum(),
            baselineIsoFlair['flag_in_FLAIR'].sum()/len(baselineIsoFlair)))
    outFile.write("\n{0} reads from baseline modified by at least one nt (internal/5'/3') in FLAIR ({1:.2%} of baseline, {2:.2%} of reads retained in FLAIR)\n".format(
            baselineIsoFlair['flag_modified_in_FLAIR'].sum(),
            baselineIsoFlair['flag_modified_in_FLAIR'].sum()/len(baselineIsoFlair),
            baselineIsoFlair['flag_modified_in_FLAIR'].sum()/baselineIsoFlair['flag_in_FLAIR'].sum()))
    outFile.write("\n{0} reads from baseline modified internally by at least one nt in FLAIR ({1:.2%} of baseline, {2:.2%} of reads retained in FLAIR, {3:.2%} of reads modified in FLAIR)\n".format(
            baselineIsoFlair['flag_modified_internal_in_FLAIR'].sum(),
            baselineIsoFlair['flag_modified_internal_in_FLAIR'].sum()/len(baselineIsoFlair),
            baselineIsoFlair['flag_modified_internal_in_FLAIR'].sum()/baselineIsoFlair['flag_in_FLAIR'].sum(),
            baselineIsoFlair['flag_modified_internal_in_FLAIR'].sum()/baselineIsoFlair['flag_modified_in_FLAIR'].sum()))
    # Number and % modified in both isoseq3 cluster and FLAIR
    outFile.write("\n{0} reads from baseline modified by at least one nt (internal/5'/3') in isoseq3 cluster and FLAIR ({1:.2%} of baseline, {2:.2%} of reads retained in isoseq3 cluster, {3:.2%} of reads retained in FLAIR)\n".format(
            baselineIsoFlair['flag_modified_in_isoseq3Cluster_FLAIR'].sum(),
            baselineIsoFlair['flag_modified_in_isoseq3Cluster_FLAIR'].sum()/len(baselineIsoFlair),
            baselineIsoFlair['flag_modified_in_isoseq3Cluster_FLAIR'].sum()/baselineIsoFlair['flag_in_isoseq3Cluster'].sum(),
            baselineIsoFlair['flag_modified_in_isoseq3Cluster_FLAIR'].sum()/baselineIsoFlair['flag_in_FLAIR'].sum()))
    outFile.write("\n{0} reads from baseline modified internally by at least one nt in isoseq3 cluster and FLAIR ({1:.2%} of baseline, {2:.2%} of reads retained in isoseq3 cluster, {3:.2%} of reads retained in FLAIR, {4:.2%} of reads modified in isoseq3 cluster, {5:.2%} of reads modified in FLAIR, {6:.2%} of reads modified internally in isoseq3 cluster, {7:.2%} of reads modified internally in FLAIR)\n".format(
            baselineIsoFlair['flag_modified_internal_in_isoseq3Cluster_FLAIR'].sum(),
            baselineIsoFlair['flag_modified_internal_in_isoseq3Cluster_FLAIR'].sum()/len(baselineIsoFlair),
            baselineIsoFlair['flag_modified_internal_in_isoseq3Cluster_FLAIR'].sum()/baselineIsoFlair['flag_in_isoseq3Cluster'].sum(),
            baselineIsoFlair['flag_modified_internal_in_isoseq3Cluster_FLAIR'].sum()/baselineIsoFlair['flag_in_FLAIR'].sum(),
            baselineIsoFlair['flag_modified_internal_in_isoseq3Cluster_FLAIR'].sum()/baselineIsoFlair['flag_modified_in_isoseq3Cluster'].sum(),
            baselineIsoFlair['flag_modified_internal_in_isoseq3Cluster_FLAIR'].sum()/baselineIsoFlair['flag_modified_in_FLAIR'].sum(),
            baselineIsoFlair['flag_modified_internal_in_isoseq3Cluster_FLAIR'].sum()/baselineIsoFlair['flag_modified_internal_in_isoseq3Cluster'].sum(),
            baselineIsoFlair['flag_modified_internal_in_isoseq3Cluster_FLAIR'].sum()/baselineIsoFlair['flag_modified_internal_in_FLAIR'].sum()))

    # Output sum of flags
    baselineIsoFlair.loc[:,'flag_in_baseline'] = 1
    cols = ['flag_in_baseline',
            'flag_in_SQANTI3filter',
            'flag_in_isoseq3Cluster',
            'flag_in_FLAIR',
            'flag_in_SQANTI3filter_isoseq3Cluster',
            'flag_in_SQANTI3filter_FLAIR',
            'flag_in_isoseq3Cluster_FLAIR',
            'flag_in_SQANTI3filter_isoseq3Cluster_FLAIR',
            'flag_in_SQANTI3filter_not_in_FLAIR',
            'flag_in_SQANTI3filter_not_in_isoseq3Cluster',
            'flag_in_FLAIR_not_in_SQANTI3filter',
            'flag_in_FLAIR_not_in_isoseq3Cluster',
            'flag_in_isoseq3Cluster_not_in_SQANTI3filter',
            'flag_in_isoseq3Cluster_not_in_FLAIR',
            'flag_modified_in_isoseq3Cluster',
            'flag_modified_internal_in_isoseq3Cluster',
            'flag_modified_in_FLAIR',
            'flag_modified_internal_in_FLAIR',
            'flag_modified_in_isoseq3Cluster_FLAIR',
            'flag_modified_internal_in_isoseq3Cluster_FLAIR']
    sumDF = pd.DataFrame(baselineIsoFlair[cols].sum()).T.rename(columns={
            'flag_in_baseline':'sum_in_baseline_read',
            'flag_in_SQANTI3filter':'sum_in_SQANTI3filter',
            'flag_in_isoseq3Cluster':'sum_in_isoseq3Cluster',
            'flag_in_FLAIR':'sum_in_FLAIR',
            'flag_in_SQANTI3filter_isoseq3Cluster':'sum_in_SQANTI3filter_isoseq3Cluster',
            'flag_in_SQANTI3filter_FLAIR':'sum_in_SQANTI3filter_FLAIR',
            'flag_in_isoseq3Cluster_FLAIR':'sum_in_isoseq3Cluster_FLAIR',
            'flag_in_SQANTI3filter_isoseq3Cluster_FLAIR':'sum_in_SQANTI3filter_isoseq3Cluster_FLAIR',
            'flag_in_SQANTI3filter_not_in_FLAIR':'sum_in_SQANTI3filter_not_in_FLAIR',
            'flag_in_SQANTI3filter_not_in_isoseq3Cluster':'sum_in_SQANTI3filter_not_in_isoseq3Cluster',
            'flag_in_FLAIR_not_in_SQANTI3filter':'sum_in_FLAIR_not_in_SQANTI3filter',
            'flag_in_FLAIR_not_in_isoseq3Cluster':'sum_in_FLAIR_not_in_isoseq3Cluster',
            'flag_in_isoseq3Cluster_not_in_SQANTI3filter':'sum_in_isoseq3Cluster_not_in_SQANTI3filter',
            'flag_in_isoseq3Cluster_not_in_FLAIR':'sum_in_isoseq3Cluster_not_in_FLAIR',
            'flag_modified_in_isoseq3Cluster':'sum_modified_in_isoseq3Cluster',
            'flag_modified_internal_in_isoseq3Cluster':'sum_modified_internal_in_isoseq3Cluster',
            'flag_modified_in_FLAIR':'sum_modified_in_FLAIR',
            'flag_modified_internal_in_FLAIR':'sum_modified_internal_in_FLAIR',
            'flag_modified_in_isoseq3Cluster_FLAIR':'sum_modified_in_isoseq3Cluster_FLAIR',
            'flag_modified_internal_in_isoseq3Cluster_FLAIR':'sum_modified_internal_in_isoseq3Cluster_FLAIR'})
    # Add other summary counts of interest
    sumDF['num_baseline_gene'] = baselineIsoFlair['gene_id'].nunique()
    sumDF['num_single_read_baseline_gene'] = baselineIsoFlair[baselineIsoFlair['flag_single_read_gene_baseline']==1]['gene_id'].nunique()
    readNumFlag = [flag for flag in baselineIsoFlair.columns if ("flag_at_least_" in flag) and ("_reads_per_gene" in flag)]
    if len(readNumFlag) > 0:
        readNumFlagName = readNumFlag[0][4:]
        sumDF['num_gene_'+readNumFlagName] = baselineIsoFlair[baselineIsoFlair[readNumFlag[0]]==1]['gene_id'].nunique()
        sumDF['num_gene_less_than_'+readNumFlagName[10:]] = baselineIsoFlair[baselineIsoFlair[readNumFlag[0]]==0]['gene_id'].nunique()
        sumDF['num_at_least_'+readNumFlagName[10:-15]+'_read_baseline_gene_in_SQANTI3filter'] = baselineIsoFlair[(baselineIsoFlair['flag_in_SQANTI3filter']==1)&(baselineIsoFlair[readNumFlag[0]]==1)]['gene_id'].nunique()
        sumDF['num_at_least_'+readNumFlagName[10:-15]+'_read_baseline_gene_in_isoseq3Cluster_no_mutligene_cluster'] = baselineIsoFlair[(baselineIsoFlair['flag_in_isoseq3Cluster']==1)&(baselineIsoFlair[readNumFlag[0]]==1)&(baselineIsoFlair['flag_multigene_isoseq3Cluster']==0)]['gene_id'].nunique()
        sumDF['num_at_least_'+readNumFlagName[10:-15]+'_read_baseline_gene_in_FLAIR_no_mutligene_cluster'] = baselineIsoFlair[(baselineIsoFlair['flag_in_FLAIR']==1)&(baselineIsoFlair[readNumFlag[0]]==1)&(baselineIsoFlair['flag_multigene_FLAIR']==0)]['gene_id'].nunique()
        sumDF['perc_at_least_'+readNumFlagName[10:-15]+'_read_baseline_gene_in_SQANTI3filter'] = "{0:.2%}".format(baselineIsoFlair[(baselineIsoFlair['flag_in_SQANTI3filter']==1)&(baselineIsoFlair[readNumFlag[0]]==1)]['gene_id'].nunique()/baselineIsoFlair[baselineIsoFlair[readNumFlag[0]]==1]['gene_id'].nunique())
        sumDF['perc_at_least_'+readNumFlagName[10:-15]+'_read_baseline_gene_in_isoseq3Cluster_no_mutligene_cluster'] = "{0:.2%}".format(baselineIsoFlair[(baselineIsoFlair['flag_in_isoseq3Cluster']==1)&(baselineIsoFlair[readNumFlag[0]]==1)&(baselineIsoFlair['flag_multigene_isoseq3Cluster']==0)]['gene_id'].nunique()/baselineIsoFlair[baselineIsoFlair[readNumFlag[0]]==1]['gene_id'].nunique())
        sumDF['perc_at_least_'+readNumFlagName[10:-15]+'_read_baseline_gene_in_FLAIR_no_mutligene_cluster'] = "{0:.2%}".format(baselineIsoFlair[(baselineIsoFlair['flag_in_FLAIR']==1)&(baselineIsoFlair[readNumFlag[0]]==1)&(baselineIsoFlair['flag_multigene_FLAIR']==0)]['gene_id'].nunique()/baselineIsoFlair[baselineIsoFlair[readNumFlag[0]]==1]['gene_id'].nunique())
    sumDF['num_single_read_baseline_gene'] = baselineIsoFlair[baselineIsoFlair['flag_single_read_gene_baseline']==1]['gene_id'].nunique()
    sumDF['num_multigene_cluster_in_isoseq3Cluster'] = baselineIsoFlair[(baselineIsoFlair['flag_multigene_isoseq3Cluster']==1)&(baselineIsoFlair['flag_in_isoseq3Cluster']==1)]['iso_cluster_id'].nunique()
    sumDF['num_multigene_cluster_in_FLAIR'] = baselineIsoFlair[(baselineIsoFlair['flag_multigene_FLAIR']==1)&(baselineIsoFlair['flag_in_FLAIR']==1)]['flair_transcript_id'].nunique()
    sumDF['num_single_gene_cluster_in_isoseq3Cluster'] = baselineIsoFlair[(baselineIsoFlair['flag_multigene_isoseq3Cluster']==0)&(baselineIsoFlair['flag_in_isoseq3Cluster']==1)]['iso_cluster_id'].nunique()
    sumDF['num_single_gene_cluster_in_FLAIR'] = baselineIsoFlair[(baselineIsoFlair['flag_multigene_FLAIR']==0)&(baselineIsoFlair['flag_in_FLAIR']==1)]['flair_transcript_id'].nunique()
    sumDF['num_multi_read_baseline_gene'] = baselineIsoFlair[baselineIsoFlair['flag_single_read_gene_baseline']==0]['gene_id'].nunique()
    sumDF['num_multi_read_baseline_gene_in_SQANTI3filter'] = baselineIsoFlair[(baselineIsoFlair['flag_in_SQANTI3filter']==1)&(baselineIsoFlair['flag_single_read_gene_baseline']==0)]['gene_id'].nunique()
    sumDF['num_multi_read_baseline_gene_in_isoseq3Cluster_no_mutligene_cluster'] = baselineIsoFlair[(baselineIsoFlair['flag_in_isoseq3Cluster']==1)&(baselineIsoFlair['flag_single_read_gene_baseline']==0)&(baselineIsoFlair['flag_multigene_isoseq3Cluster']==0)]['gene_id'].nunique()
    sumDF['num_multi_read_baseline_gene_in_FLAIR_no_mutligene_cluster'] = baselineIsoFlair[(baselineIsoFlair['flag_in_FLAIR']==1)&(baselineIsoFlair['flag_single_read_gene_baseline']==0)&(baselineIsoFlair['flag_multigene_FLAIR']==0)]['gene_id'].nunique()
    sumDF['perc_multi_read_baseline_gene_in_SQANTI3filter'] = "{0:.2%}".format(baselineIsoFlair[(baselineIsoFlair['flag_in_SQANTI3filter']==1)&(baselineIsoFlair['flag_single_read_gene_baseline']==0)]['gene_id'].nunique()/baselineIsoFlair[baselineIsoFlair['flag_single_read_gene_baseline']==0]['gene_id'].nunique())
    sumDF['perc_multi_read_baseline_gene_in_isoseq3Cluster_no_mutligene_cluster'] = "{0:.2%}".format(baselineIsoFlair[(baselineIsoFlair['flag_in_isoseq3Cluster']==1)&(baselineIsoFlair['flag_single_read_gene_baseline']==0)&(baselineIsoFlair['flag_multigene_isoseq3Cluster']==0)]['gene_id'].nunique()/baselineIsoFlair[baselineIsoFlair['flag_single_read_gene_baseline']==0]['gene_id'].nunique())
    sumDF['perc_multi_read_baseline_gene_in_FLAIR_no_mutligene_cluster'] = "{0:.2%}".format(baselineIsoFlair[(baselineIsoFlair['flag_in_FLAIR']==1)&(baselineIsoFlair['flag_single_read_gene_baseline']==0)&(baselineIsoFlair['flag_multigene_FLAIR']==0)]['gene_id'].nunique()/baselineIsoFlair[baselineIsoFlair['flag_single_read_gene_baseline']==0]['gene_id'].nunique())
    
    print(sumDF.to_csv(index=False, header=False),end="")


if __name__ == '__main__':
    # Parse command line arguments
    global args
    args = getOptions()
    main()
