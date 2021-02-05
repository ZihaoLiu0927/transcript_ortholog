#!/usr/bin/env python

import argparse
import pandas as pd
import numpy as np
import sqlite3
import csv

def getOptions():
    # Parse command line arguments
    parser = argparse.ArgumentParser(description="")

    # Input data
    parser.add_argument("-b", "--baseline", dest="inBase", required=True, help="Baseline GTF file")
    parser.add_argument("-i", "--isoseq3Cluster", dest="inIso", required=True, help="Isoseq3 clustered, re-aligned GTF file")
    parser.add_argument("-r", "--report", dest="inReport", required=True, help="Isoseq3 clustered report to match baseline read ID values with isoseq3 cluster ID values")
    parser.add_argument("-s", "--SQANTIfilter", dest="inSQANTI", required=True, help="SQANTI filter GTF file")
    parser.add_argument("-f", "--FLAIR", dest="inFLAIR", required=True, help="FLAIR GTF file")
    parser.add_argument("-m", "--FLAIR-map", dest="inMap", required=True, help="FLAIR read map file")
    parser.add_argument("-n", "--number-reads", dest="inNum", default=10, type=int, choices=range(0,1000000), help="Number of reads to be used for flagging reads in genes with at least N reads evidence and output associated gene values for downstream analysis")

    # Output data
    parser.add_argument("-o", "--output-file", dest="outFile", required=True, help="Path to output CSV file of flags")
    parser.add_argument("-d", "--output-dir", dest="outDir", required=False, help="Directory for modified read GTF files")
    parser.add_argument("-g", "--gene-list", dest="outGene", required=False, help="Path to output file of genes with at least N reads (set by -n argument) of evidence")
    parser.add_argument("-db", "--db-file", dest="tempDB", required=False, help="Temporary SQL database file to reduce memory requiredfor SQL merge", default=":memory:")

    args = parser.parse_args()
    return args

def get_internal_exonID(baseXcrpt):
    baseXcrpt['first_exonID'] = baseXcrpt['exonID'].str.split("|").str[0]
    baseXcrpt['last_exonID'] = baseXcrpt['exonID'].str.split("|").str[-1]
    baseXcrpt['internal_exonID'] = np.where(baseXcrpt['first_exonID']==baseXcrpt['last_exonID'], "monoexon",
             np.where(baseXcrpt['exonID'].str.split("|").str[1:-1].str.join("|")=="",
             baseXcrpt['first_exonID'].str.split("_").str[-1] + "|" + baseXcrpt['last_exonID'].str.split("_").str[:-1].str.join("_"),
             np.where(baseXcrpt['exonID'].str.split("|").str[1:-1].str.join("|")!="",
                      baseXcrpt['first_exonID'].str.split("_").str[-1] + "|" + baseXcrpt['exonID'].str.split("|").str[1:-1].str.join("|") + "|" + baseXcrpt['last_exonID'].str.split("_").str[:-1].str.join("_"),"oops"))) 
    return baseXcrpt.drop(columns=['first_exonID','last_exonID'])


def main():
    # Get all input files
    baseGTF = pd.read_csv(args.inBase,names=['chr','source','feature','start','end','score',
                                        'strand','frame','attribute'], sep="\t",low_memory=False)
    isoGTF = pd.read_csv(args.inIso,names=['chr','source','feature','start','end','score',
                                        'strand','frame','attribute'], sep="\t",low_memory=False)
    isoReport = pd.read_csv(args.inReport, low_memory=False)
    sqantiGTF = pd.read_csv(args.inSQANTI,names=['chr','source','feature','start','end','score',
                                        'strand','frame','attribute'], sep="\t",low_memory=False)
    flairGTF = pd.read_csv(args.inFLAIR,names=['chr','source','feature','start','end','score',
                                        'strand','frame','attribute'], sep="\t",low_memory=False)
    flairMap = pd.read_csv(args.inMap,names=['isoform','readID'],sep="\t",low_memory=False)
    
    # For all GTF files, retain only exon features and make exonID (exon_chr_start_end)
    baseGTF = baseGTF[baseGTF['feature']=="exon"]
    baseGTF['exonID'] = baseGTF['feature']+"_"+baseGTF['chr'].map(str)+"_"+baseGTF['start'].map(str)+"_"+baseGTF['end'].map(str)
    isoGTF = isoGTF[isoGTF['feature']=="exon"]
    isoGTF['exonID'] = isoGTF['feature']+"_"+isoGTF['chr'].map(str)+"_"+isoGTF['start'].map(str)+"_"+isoGTF['end'].map(str)
    sqantiGTF = sqantiGTF[sqantiGTF['feature']=="exon"]
    sqantiGTF['exonID'] = sqantiGTF['feature']+"_"+sqantiGTF['chr'].map(str)+"_"+sqantiGTF['start'].map(str)+"_"+sqantiGTF['end'].map(str)
    flairGTF = flairGTF[flairGTF['feature']=="exon"]
    flairGTF['exonID'] = flairGTF['feature']+"_"+flairGTF['chr'].map(str)+"_"+flairGTF['start'].map(str)+"_"+flairGTF['end'].map(str)

    # Get readID values from baseline and SQANTI3 filter GTF transcript_id values
    baseGTF['readID'] = baseGTF['attribute'].str.split(" ").str[1].str.strip("\";")
    sqantiGTF['readID'] = sqantiGTF['attribute'].str.split(" ").str[3].str.strip("\";")

    # Get cluster_id values from isoseq3 cluster GTF transcript_id values
    isoGTF['cluster_id'] = isoGTF['attribute'].str.split(" ").str[1].str.strip("\";")
    
    # Get transcript_id values from FLAIR GTF (could be readID or isoform)
    flairGTF['transcript_id'] = flairGTF['attribute'].str.split(" ").str[3].str.strip("\";")
    
    # Get gene_id values from GTF
    baseGTF['gene_id'] = baseGTF['attribute'].str.split(" ").str[3].str.strip("\";")
    sqantiGTF['gene_id'] = sqantiGTF['attribute'].str.split(" ").str[1].str.strip("\";")
    isoGTF['gene_id'] = isoGTF['attribute'].str.split(" ").str[3].str.strip("\";")
    flairGTF['gene_id'] = flairGTF['attribute'].str.split(" ").str[1].str.strip("\";")
    
    # Groupby readID/cluster_id and get concatenated list of exonID values
    baseXcrpt = baseGTF.sort_values(by=['start'],ascending=True).groupby('readID').agg({
            'chr':'first',
            'start':'min',
            'end':'max',
            'strand':'first',
            'gene_id':'first',
            'exonID':lambda x: "|".join(x)}).reset_index()
    sqantiXcrpt = sqantiGTF.sort_values(by=['start'],ascending=True).groupby('readID').agg({
            'chr':'first',
            'start':'min',
            'end':'max',
            'strand':'first',
            'gene_id':'first',
            'exonID':lambda x: "|".join(x)}).reset_index()
    isoXcrpt = isoGTF.sort_values(by=['start'],ascending=True).groupby('cluster_id').agg({
            'chr':'first',
            'start':'min',
            'end':'max',
            'strand':'first',
            'gene_id':'first',
            'exonID':lambda x: "|".join(x)}).reset_index()
    flairXcrpt = flairGTF.sort_values(by=['start'],ascending=True).groupby('transcript_id').agg({
            'chr':'first',
            'start':'min',
            'end':'max',
            'strand':'first',
            'gene_id':'first',
            'exonID':lambda x: "|".join(x)}).reset_index()
    # Get internal splice sites/exons
    baseXcrpt = get_internal_exonID(baseXcrpt)
#    sqantiXcrpt = get_internal_exonID(sqantiXcrpt)
    isoXcrpt = get_internal_exonID(isoXcrpt)
    flairXcrpt = get_internal_exonID(flairXcrpt)
   
    # Flag if baseline readID is present in and SQANTI3 readID list
    # (do not need to compare coordinates since SQANTI3 filter does not change)
    baseXcrpt['flag_in_SQANTI3filter'] = np.where(baseXcrpt['readID'].isin(sqantiXcrpt['readID']),1,0)

    # Merge isoseq3 cluster GTF derived cluster_id values with the cluster report
    #     by cluster_id to get the readID values associated with each cluster_id
    # Output will contain NA values in reads that are not represented in cluster_id
    con = sqlite3.connect(args.tempDB)
    cur = con.cursor()
    isoXcrpt.to_sql("isoXcrpt", con, if_exists="replace")
    isoReport.to_sql("report", con, if_exists="replace")
    cur.execute("CREATE TABLE mergeIso AS SELECT in1.read_id, in1.read_type, in1.cluster_id, in2.* "
                "FROM report in1 LEFT JOIN isoXcrpt in2 "
                "ON in1.cluster_id = in2.cluster_id")
    isoDF = pd.read_sql("SELECT * FROM mergeIso", con).drop(columns=['index','cluster_id:1']).rename(
            columns={'read_id':'readID'})
    con.close()
    # Add prefix to isoseq3 cluster columns
    isoDF.columns = "iso_"+isoDF.columns
    # Merge baseline and isoseq3 cluster
    con = sqlite3.connect(args.tempDB)
    cur = con.cursor()
    isoDF.to_sql("iso", con, if_exists="replace")
    baseXcrpt.to_sql("baseline", con, if_exists="replace")
    cur.execute("CREATE TABLE baseIso AS SELECT in1.*, in2.* "
                "FROM baseline in1 LEFT JOIN iso in2 "
                "ON in1.readID = in2.iso_readID")
    baselineIso = pd.read_sql("SELECT * from baseIso", con).drop(columns=['index','index:1'])
    con.close()
    
    # Flag reads that were retained in isoseq3 cluster
    baselineIso['flag_in_isoseq3Cluster'] = np.where(~baselineIso['iso_exonID'].isna(),1,0)
    
    # Flag reads that were modified in isoseq3 cluster (reads not retained are 0)
    # Any modification (internal/5'/3')
    baselineIso['flag_modified_in_isoseq3Cluster'] = np.where((~baselineIso['iso_exonID'].isna())&
               (baselineIso['exonID']!=baselineIso['iso_exonID']),1,0)
    # Flag reads that were modified internally in isoseq3 cluster (reads not retained are 0)
    # Only internal modification (no 5'/3' change)
    baselineIso['flag_modified_internal_in_isoseq3Cluster'] = np.where((~baselineIso['iso_internal_exonID'].isna())&
               (baselineIso['internal_exonID']!=baselineIso['iso_internal_exonID']),1,0)
    
    # Flag reads that were retained in both isoseq3 cluster and SQANTI3 filter
    baselineIso['flag_in_SQANTI3filter_isoseq3Cluster'] = np.where(
            (baselineIso['flag_in_SQANTI3filter']==1)&(baselineIso['flag_in_isoseq3Cluster']==1),
            1,0)
    
    # Split FLAIR readID values by "," to get stacked instead of list
    splitList = flairMap['readID'].str.split(",").apply(pd.Series, 1).stack()
    splitList.index = splitList.index.droplevel(-1)
    tempDF = flairMap.copy()
    del(tempDF['readID'])
    splitMap = tempDF.join(splitList.rename('readID'))
    # Split FLAIR isoform names into transcript_id and gene_id
    splitMap['transcript_id'] = splitMap['isoform'].str.rsplit("_",1).str[0]
    splitMap['gene_id'] = splitMap['isoform'].str.rsplit("_",1).str[1]
    
    # Merge FLAIR GTF with file of read ID values associated with isoforms
    con = sqlite3.connect(args.tempDB)
    cur = con.cursor()
    flairXcrpt.to_sql("flairXcrpt", con, if_exists="replace")
    splitMap.to_sql("splitMap", con, if_exists="replace")
    cur.execute("CREATE TABLE mergeFlair AS SELECT in1.transcript_id, in1.readID, in2.* "
                "FROM splitMap in1 LEFT JOIN flairXcrpt in2 "
                "ON in1.transcript_id = in2.transcript_id")
    flairDF = pd.read_sql("SELECT * from mergeFlair", con).drop(columns=['index','transcript_id:1'])
    con.close()
    
    # Add prefix to FLAIR columns
    flairDF.columns = "flair_"+flairDF.columns
    
    # Merge baseline with isoseeq3 cluster and FLAIR
    con = sqlite3.connect(args.tempDB)
    cur = con.cursor()
    flairDF.to_sql("flair", con, if_exists="replace")
    baselineIso.to_sql("baselineIso", con, if_exists="replace")
    cur.execute("CREATE TABLE baseIsoFLAIR AS SELECT in1.*, in2.* "
                "FROM baselineIso in1 LEFT JOIN flair in2 "
                "ON in1.readID = in2.flair_readID")
    baselineIsoFlair = pd.read_sql("SELECT * from baseIsoFLAIR", con).drop(columns=['index','index:1'])
    con.close()
    
    # Flag reads that were retained in FLAIR
    baselineIsoFlair['flag_in_FLAIR'] = np.where(~baselineIsoFlair['flair_exonID'].isna(),1,0)
    # Flag reads that were modified in FLAIR (reads not retained are 0)
    # Any modification (internal/5'/3')
    baselineIsoFlair['flag_modified_in_FLAIR'] = np.where((~baselineIsoFlair['flair_exonID'].isna())&
               (baselineIsoFlair['exonID']!=baselineIsoFlair['flair_exonID']),1,0)
    # Flag reads that were internally modified in FLAIR (reads not retained are 0)
    # Only internal modification (no 5'/3' change)
    baselineIsoFlair['flag_modified_internal_in_FLAIR'] = np.where((~baselineIsoFlair['flair_internal_exonID'].isna())&
               (baselineIsoFlair['internal_exonID']!=baselineIsoFlair['flair_internal_exonID']),1,0)
    # Flag reads that were retained in both SQANTI3 filter and FLAIR
    baselineIsoFlair['flag_in_SQANTI3filter_FLAIR'] = np.where(
            (baselineIsoFlair['flag_in_SQANTI3filter']==1)&(baselineIsoFlair['flag_in_FLAIR']==1),
            1,0)
    # Flag reads that were retained in both isoseq3 cluster and FLAIR
    baselineIsoFlair['flag_in_isoseq3Cluster_FLAIR'] = np.where(
            (baselineIsoFlair['flag_in_isoseq3Cluster']==1)&(baselineIsoFlair['flag_in_FLAIR']==1),
            1,0)
    # Flag reads that were retained in SQANTI3 filter, isoseq3 cluster and FLAIR
    baselineIsoFlair['flag_in_SQANTI3filter_isoseq3Cluster_FLAIR'] = np.where(
            (baselineIsoFlair['flag_in_SQANTI3filter_isoseq3Cluster']==1)&(baselineIsoFlair['flag_in_FLAIR']==1),
            1,0)
    # Flag reads that were modified in isoseq3 cluster and FLAIR (reads not retained are 0)
    # Any modification (internal/5'/3')
    baselineIsoFlair['flag_modified_in_isoseq3Cluster_FLAIR'] = np.where(
            (baselineIsoFlair['flag_modified_in_isoseq3Cluster']==1)&(baselineIsoFlair['flag_modified_in_FLAIR']==1),
            1,0)
    # Flag reads that were internally modified in isoseq3 cluster FLAIR (reads not retained are 0)
    # Only internal modification (no 5'/3' change)
    baselineIsoFlair['flag_modified_internal_in_isoseq3Cluster_FLAIR'] = np.where(
            (baselineIsoFlair['flag_modified_internal_in_isoseq3Cluster']==1)&(baselineIsoFlair['flag_modified_internal_in_FLAIR']==1),
            1,0)
    # Flag reads that were retained in SQANTI3 filter but not in FLAIR
    baselineIsoFlair['flag_in_SQANTI3filter_not_in_FLAIR'] = np.where(
            (baselineIsoFlair['flag_in_SQANTI3filter']==1)&(baselineIsoFlair['flag_in_FLAIR']==0),
            1,0)
    # Flag reads that were retained in SQANTI3 filter but not in isoseq3 cluster
    baselineIsoFlair['flag_in_SQANTI3filter_not_in_isoseq3Cluster'] = np.where(
            (baselineIsoFlair['flag_in_SQANTI3filter']==1)&(baselineIsoFlair['flag_in_isoseq3Cluster']==0),
            1,0)
    # Flag reads that were retained in FLAIR but not in SQANTI3 filter
    baselineIsoFlair['flag_in_FLAIR_not_in_SQANTI3filter'] = np.where(
            (baselineIsoFlair['flag_in_SQANTI3filter']==0)&(baselineIsoFlair['flag_in_FLAIR']==1),
            1,0)
    # Flag reads that were retained in FLAIR but not in isoseq3 cluster
    baselineIsoFlair['flag_in_FLAIR_not_in_isoseq3Cluster'] = np.where(
            (baselineIsoFlair['flag_in_isoseq3Cluster']==0)&(baselineIsoFlair['flag_in_FLAIR']==1),
            1,0)
    # Flag reads that were retained in isoseq3 cluster but not in SQANTI3 filter
    baselineIsoFlair['flag_in_isoseq3Cluster_not_in_SQANTI3filter'] = np.where(
            (baselineIsoFlair['flag_in_isoseq3Cluster']==1)&(baselineIsoFlair['flag_in_SQANTI3filter']==0),
            1,0)
    # Flag reads that were retained in isoseq3 cluster but not in FLAIR
    baselineIsoFlair['flag_in_isoseq3Cluster_not_in_FLAIR'] = np.where(
            (baselineIsoFlair['flag_in_isoseq3Cluster']==1)&(baselineIsoFlair['flag_in_FLAIR']==0),
            1,0)
    # Flag reads/genes that are a single read in a gene
    baselineIsoFlair['baseline_read_per_gene'] = baselineIsoFlair.groupby('gene_id')['readID'].transform('count')
    baselineIsoFlair['flag_single_read_gene_baseline'] = np.where(baselineIsoFlair['baseline_read_per_gene']==1,1,0)
    # Flag reads associated with genes that have at least N reads of evidence 
    baselineIsoFlair['flag_at_least_'+str(args.inNum)+'_reads_per_gene'] = np.where(
            baselineIsoFlair['baseline_read_per_gene']>=args.inNum,1,0)
    if args.outGene is not None:
        # Output list of genes with at least N reads of evidence
        baselineIsoFlair[baselineIsoFlair['flag_at_least_'+str(args.inNum)+'_reads_per_gene']==1]['gene_id'].drop_duplicates().to_csv(args.outGene, index=False, header=False)
    # Count number of reads per cluster in isoseq3 cluster and FLAIR and flag single read clusters
    baselineIsoFlair['iso_readID_per_cluster'] = baselineIsoFlair.groupby('iso_cluster_id')['readID'].transform('count')
    baselineIsoFlair['flair_readID_per_cluster'] = baselineIsoFlair.groupby('flair_transcript_id')['readID'].transform('count')
    baselineIsoFlair['flag_single_read_cluster_isoseq3Cluster'] = np.where(baselineIsoFlair['iso_readID_per_cluster']==1,1,0)
    baselineIsoFlair['flag_single_read_cluster_FLAIR'] = np.where(baselineIsoFlair['flair_readID_per_cluster']==1,1,0)
    # Flag monoexon reads/clusters
    baselineIsoFlair['flag_monoexon_baseline'] = np.where(baselineIsoFlair['internal_exonID']=="monoexon",1,0)
    baselineIsoFlair['flag_monoexon_isoseq3Cluster'] = np.where(baselineIsoFlair['iso_internal_exonID']=="monoexon",1,0)
    baselineIsoFlair['flag_monoexon_FLAIR'] = np.where(baselineIsoFlair['flair_internal_exonID']=="monoexon",1,0)
    # Flag clusters in isoseq3 cluster or FLAIR that contain baseline reads from multiple genes
    baselineIsoFlair['iso_num_gene_in_cluster'] = baselineIsoFlair.groupby('iso_cluster_id')['gene_id'].transform('nunique')
    baselineIsoFlair['flair_num_gene_in_cluster'] = baselineIsoFlair.groupby('flair_transcript_id')['gene_id'].transform('nunique')
    baselineIsoFlair['flag_multigene_isoseq3Cluster'] = np.where(baselineIsoFlair['iso_num_gene_in_cluster']>1,1,0)
    baselineIsoFlair['flag_multigene_FLAIR'] = np.where(baselineIsoFlair['flair_num_gene_in_cluster']>1,1,0)
    # Flag genes that are single cluster genes in isoseq3 cluster or FLAIR (excluding multigene clusters)
#    baselineIsoFlair['iso_single_gene_cluster_per_gene'] = baselineIsoFlair.groupby('gene_id')['iso_cluster_id'].transform('nunique')
#    baselineIsoFlair['flair_single_gene_cluster_per_gene'] = baselineIsoFlair[(baselineIsoFlair['flag_multigene_FLAIR']==0)&(baselineIsoFlair['flag_in_FLAIR']==1)].groupby('gene_id')['flair_transcript_id'].transform('nunique')
#    baselineIsoFlair['flag_one_single_gene_cluster_per_gene_isoseq3Cluster'] = np.where(baselineIsoFlair['iso_single_gene_cluster_per_gene']==1,1,0)
#    baselineIsoFlair['flag_one_single_gene_cluster_per_gene_FLAIR'] = np.where(baselineIsoFlair['flair_single_gene_cluster_per_gene']==1,1,0)

    # Output flag file
    baselineIsoFlair.to_csv(args.outFile, index=False)

    # If GTF output directory given
    # Output GTF files for each read that is modified in isoseq3 cluster and/or FLAIR
    if args.outDir is not None:
        for index,row in baselineIsoFlair[(baselineIsoFlair['flag_modified_internal_in_isoseq3Cluster']==1)|(baselineIsoFlair['flag_modified_internal_in_FLAIR']==1)].iterrows():
            # Open new file for readID
            #     (replace any "/" in readID iwth "_")
            tempOut = open(args.outDir+"/"+row['readID'].replace("/","_")+".gtf",'w')
            tempBase = baseGTF[baseGTF['readID']==row['readID']].reset_index(drop=True)
            tempBase['readID'] = "baseline_"+tempBase['readID']
            tempBase['new_attribute'] = "transcript_id \""+tempBase['readID']+"\"; gene_id \""+tempBase['gene_id']+"\";"
            tempOut.write(tempBase[['chr','source','feature','start','end','score','strand','frame',
                                   'new_attribute']].to_csv(sep="\t",index=False,header=False,
                                    doublequote=False,quoting=csv.QUOTE_NONE))
            if row['flag_modified_in_isoseq3Cluster'] == 1:
                tempIso = isoGTF[isoGTF['cluster_id']==row['iso_cluster_id']].reset_index(drop=True)
                tempIso['cluster_id'] = "iso_"+tempIso['cluster_id']
                tempIso['new_attribute'] = "transcript_id \""+tempIso['cluster_id']+"\"; gene_id \""+tempIso['gene_id']+"\";"
                tempOut.write(tempIso[['chr','source','feature','start','end','score','strand','frame',
                                       'new_attribute']].to_csv(sep="\t",index=False,header=False,
                                        doublequote=False,quoting=csv.QUOTE_NONE))
            if row['flag_modified_in_FLAIR'] == 1:
                tempFLAIR = flairGTF[flairGTF['transcript_id']==row['flair_transcript_id']].reset_index(drop=True)
                tempFLAIR['transcript_id'] = "flair_"+tempFLAIR['transcript_id']
                tempFLAIR['new_attribute'] = "transcript_id \""+tempFLAIR['transcript_id']+"\"; gene_id \""+tempFLAIR['gene_id']+"\";"
                tempOut.write(tempFLAIR[['chr','source','feature','start','end','score','strand','frame',
                                         'new_attribute']].to_csv(sep="\t",index=False,header=False,
                                          doublequote=False,quoting=csv.QUOTE_NONE))
            tempOut.close()
        
if __name__ == '__main__':
    # Parse command line arguments
    global args
    args = getOptions()
    main()
