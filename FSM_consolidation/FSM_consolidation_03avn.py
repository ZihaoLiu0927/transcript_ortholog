#!/usr/bin/env python

import argparse
import pandas as pd
import numpy as np
import sqlite3

# Import custon functions for FSM consolidation
import FSM_consolidation_functions as FSM

def getOptions():
    # Parse command line arguments
    parser = argparse.ArgumentParser(description="Consolidate transcriptome by selecting longest (in nt) FSM and consolidating mono-exon genes to a single transcript")

    # Input data
    parser.add_argument("-i", "--input-transcripts", dest="inXcrptGene", required=True, help="CSV file with unique pairs of transcript_id to gene_id")
    parser.add_argument("-f", "--fragments", dest="inFrag", required=True, help="CSV file of fragments annotations including transcript_id column (separated by '|' if there are multiple (*_exon_fragment_annotations.csv)")
    parser.add_argument("-j", "--junctions", dest="inJunc", required=True, help="CSV file of junction annotations (*_annotated_junctions.csv)")
    parser.add_argument("-t", "--transcript-prefix", dest="inXcrptPrefix", required=True, default="tr", help="Prefix for representative transcript_id values in output (default: 'tr')")

    # Output data
    parser.add_argument("-d", "--output-directory", dest="outDir", required=True, help="Output directory")
    parser.add_argument("-p", "--output-prefix", dest="outPrefix", required=True, help="Prefix for output files")

    args = parser.parse_args()
    return args

def main():
    # Get input file of transcript_id to gene_id
    #     (will be used to get gene_id for monoexon transcripts)
    xcrptGene = pd.read_csv(args.inXcrptGene, low_memory=False)
    # Get input file of fragments
    fragDF = pd.read_csv(args.inFrag, low_memory=False)
    # Add fragment length
    fragDF['fragment_length'] = fragDF['fragment_stop'] - fragDF['fragment_start']
    
    # Get input file of junctions (will only include spliced transcripts)
    #      junction_id is not unique (pair of junction_id,transcript_id is unique)
    juncDF = pd.read_csv(args.inJunc, low_memory=False)
    # Drop empty junction_id values if present
    juncDF = juncDF.dropna()
    # Drop exon ID junction_id and set junction_coordinates as junction_id
    #     (matches the junction_id values in *_junctions_full_annotation.csv)
    juncDF = juncDF.drop(columns=['junction_id'])
    juncDF = juncDF.rename(columns={'junction_coordinates':'junction_id'})
    
    # Get donor_stop and acceptor_start from coordinates in junction_id
    juncDF['donor_stop'] = juncDF['junction_id'].str.split(":").str[1].astype(int)
    juncDF['acceptor_start'] = juncDF['junction_id'].str.split(":").str[2].astype(int)
    
    # Print out counts from input files
    print("{} transcripts and {} genes in input files".format(xcrptGene['transcript_id'].nunique(),xcrptGene['gene_id'].nunique()))
    print("{} unique junctions and {} unique fragments present in annotation".format(juncDF['junction_id'].nunique(),fragDF['fragment_id'].nunique()))
        
    # Split transcript_id in by pipes and keep all other values the same
    # Some events are found in multiple transcripts (common or constitutive) and
    #     have a piped list of transcript_id values in the transcript_id column
    splitFragDF = FSM.split_transcript_id(df=fragDF,sort_list=['transcript_id','fragment_start'])
    del(fragDF)

    # Make transcript-level junction variables by collapsing by transcript_id
    # Sort junctions by transcript_id and donor_stop
    juncDF = juncDF.sort_values(['transcript_id','donor_stop'])
    collapseJuncDF = FSM.collapse_junctions(juncDF)
    del(juncDF)
    
    # Merge transcript-level junction (multiexon transcripts only) and
    #     gene_id values (all transcripts including multi-exon) variables to get
    #     mono-exon transcripts
    con = sqlite3.connect(':memory:')
    cur = con.cursor()
    xcrptGene.to_sql("gene", con, if_exists="replace")
    collapseJuncDF.to_sql("junction", con, if_exists="replace")
    # Merge transcript level junction variables and gene_id by transcript_id
    cur.execute("CREATE TABLE merge1 AS SELECT in1.gene_id, in1.transcript_id, in2.junctionID_order "
                "FROM gene in1 LEFT JOIN junction in2 "
                "ON in1.transcript_id = in2.transcript_id")
    collapseJuncXcrptDF = pd.read_sql("SELECT * FROM merge1", con)
    # Verify merge is good
    if len(collapseJuncXcrptDF) != len(xcrptGene):
        print("WARNING: Unexpected merge")
    if len(collapseJuncXcrptDF) != len(collapseJuncDF) + len(collapseJuncXcrptDF[collapseJuncXcrptDF['junctionID_order'].isna()]):
        print("WARNING: Unexpected merge")
    del(collapseJuncDF)
    # Merge fragment-level variables and gene_id values variables to get proper gene_id value
    #     in multiexon fragments
    splitFragDF.to_sql("fragment", con, if_exists="replace")
    cur.execute("CREATE TABLE merge2 AS SELECT in1.transcript_id,in1.fragment_id,"
                "in1.chr,in1.fragment_start,in1.fragment_stop,in1.fragment_length, in2.gene_id "
                "FROM fragment in1 LEFT JOIN gene in2 "
                "ON in1.transcript_id = in2.transcript_id")
    xcrptFragDF = pd.read_sql("SELECT * FROM merge2", con)
    # Verify merge is good
    if len(xcrptFragDF) != len(splitFragDF):
        print("WARNING: Unexpected merge")
    if len(xcrptFragDF[xcrptFragDF['gene_id'].isna()]) != 0:
        print("WARNING: Unexpected merge")
    con.close()
    del(splitFragDF)

    # Flag monoexon transcripts (no junctionID_order present)
    collapseJuncXcrptDF['flag_monoexon_transcript'] = np.where(collapseJuncXcrptDF['junctionID_order'].isna(),1,0)

    # Transcripts with no junctionID_order will be set to transcript_id value
    collapseJuncXcrptDF = collapseJuncXcrptDF.fillna(-1)
    collapseJuncXcrptDF.loc[collapseJuncXcrptDF['junctionID_order']==-1,'junctionID_order'] = collapseJuncXcrptDF['transcript_id']

    # Make groups of transcripts with same junctionID_order within the same gene
    # Group names will be piped list of transcript IDs that share the same junctions
    collapseJuncXcrptDF['transcriptID_cat'] = collapseJuncXcrptDF.groupby(['junctionID_order','gene_id'])['transcript_id'].transform(func=lambda x: "|".join(x))

    # Split monoexon genes from multiexon genes check if number of transcripts
    #     in gene matches number of monoexon transcripts in gene
    collapseJuncXcrptDF['transcript_per_gene'] = collapseJuncXcrptDF.groupby('gene_id')['transcript_id'].transform('count')
    collapseJuncXcrptDF['monoexon_transcript_per_gene'] = collapseJuncXcrptDF.groupby('gene_id')['flag_monoexon_transcript'].transform('sum')
    collapseJuncXcrptDF['flag_monoexon_gene'] = np.where(collapseJuncXcrptDF['transcript_per_gene']==collapseJuncXcrptDF['monoexon_transcript_per_gene'],1,0)
    collapseMonoexonJunc = collapseJuncXcrptDF.loc[collapseJuncXcrptDF['flag_monoexon_gene']==1]
    collapseMultiexonJunc = collapseJuncXcrptDF.loc[collapseJuncXcrptDF['flag_monoexon_gene']==0]
    if len(collapseJuncXcrptDF) != len(collapseMonoexonJunc)+len(collapseMultiexonJunc):
        print("WARNING: Incorrect splitting of mono-exon and multi-exon genes")
    del(collapseJuncXcrptDF)
    
    # Split fragment dataframe by transcripts in mono-exon genes and multi-exon genes
    monoexonFrag = xcrptFragDF.loc[xcrptFragDF['transcript_id'].isin(collapseMonoexonJunc['transcript_id'])]
    multiexonFrag = xcrptFragDF.loc[xcrptFragDF['transcript_id'].isin(collapseMultiexonJunc['transcript_id'])]
    if len(xcrptFragDF) != len(monoexonFrag)+len(multiexonFrag):
        print("WARNING: Incorrect splitting of fragments from mono-exon and multi-exon genes")
    del(xcrptFragDF)
    
    # Collapse fragments of mono-exon genes by gene_id to create a single
    #     transcripts to represent each gene
    # Set transcriptID_cat and junctionID_order as piped list of all monoexon
    #     transcript_id represented in each gene
    collapseMonoexon = FSM.collapse_fragments(df=monoexonFrag,sort_list=['fragment_start'],level="gene")
    collapseMonoexon['junctionID_order'] = collapseMonoexon['transcriptID_cat']
    if len(collapseMonoexon) != collapseMonoexonJunc['gene_id'].nunique():
        print("WARNING: Incorrect collapsing of mono-exon genes")
    del(monoexonFrag,collapseMonoexonJunc)

    # Make transcript-level fragment variables by collapsing by transcript_id
    #     for multi-exon genes
    collapseMultiexonFrag = FSM.collapse_fragments(df=multiexonFrag,sort_list=['transcript_id','fragment_start'],level="transcript")
    del(multiexonFrag)
    
    # Merge transcript-level multi-exon junction and fragment variables by transcript_id
    # Verify that transcript_id is unique in both dataframes to be merged
    if len(collapseMultiexonJunc) != collapseMultiexonJunc['transcript_id'].nunique():
        print("WARNING: transcript_id not unique in collapsed junction dataframe")
    if len(collapseMultiexonFrag) != collapseMultiexonFrag['transcript_id'].nunique():
        print("WARNING: transcript_id not unique in collapsed transcript dataframe")
    con = sqlite3.connect(':memory:')
    cur = con.cursor()
    collapseMultiexonFrag.to_sql("fragment", con, if_exists="replace")
    collapseMultiexonJunc.to_sql("junction", con, if_exists="replace")
    cur.execute("CREATE TABLE merge AS SELECT in1.*, in2.* "
                "FROM fragment in1 LEFT JOIN junction in2 "
                "ON in1.transcript_id = in2.transcript_id")
    mergeMultiexon = pd.read_sql("SELECT * FROM merge", con).drop(columns=[
            'index','index:1','transcript_id:1','gene_id:1'])
    con.close()
    # Verify merge is good
    if len(collapseMultiexonFrag) != len(mergeMultiexon) or len(collapseMultiexonJunc) != len(mergeMultiexon):
        print("WARNING: Unexpected merge")
    del(collapseMultiexonFrag,collapseMultiexonJunc)
    
    # Prepare mono-exon and multi-exon dataframes for concatenation:
    #     In multi-exon drop flag_monoexon_transcript, flag_monoexon_gene,
    #         transcript_per_gene, monoexon_transcript_per_gene
    #     In mono-exon set transcript_id to transcriptID_cat value
    mergeMultiexon = mergeMultiexon.drop(
            columns=['flag_monoexon_transcript','flag_monoexon_gene',
                     'transcript_per_gene','monoexon_transcript_per_gene'])
    collapseMonoexon['transcript_id'] = collapseMonoexon['transcriptID_cat']
    
    # Concatenate mono-exon and multi-exon genes
    mergeAllDF = pd.concat([collapseMonoexon,mergeMultiexon],sort=True).reset_index(drop=True)
    del(collapseMonoexon,mergeMultiexon)

    # Get minimum start and maximum end for longest 5' and 3' ends within each group
    #     with the same junctionID_order and gene_id
    mergeAllDF['min_transcript_start'] = mergeAllDF.groupby(['junctionID_order','gene_id'])['start'].transform('min')
    mergeAllDF['max_transcript_end'] = mergeAllDF.groupby(['junctionID_order','gene_id'])['end'].transform('max')
 
    # Select the longest transcript for each junctionID_order within the same gene
    longestDF = mergeAllDF.loc[mergeAllDF.groupby(['junctionID_order','gene_id'])['transcript_length'].idxmax()]

    # Generate new transcript_id for each transcript:
    #     [prefix]_[gene_id]_# where #={1,...,n} for all n transcripts in the gene
    # First sort by gene_id and transcript_length
    longestDF = longestDF.sort_values(by=['gene_id','transcript_length'],ascending=[True,False])
    longestDF['transcript_rank_in_gene'] = longestDF.groupby('gene_id'
             )['transcript_length'].rank(method='first')
    longestDF['FSM_consolidation_transcript_id'] = args.inXcrptPrefix+"_"+longestDF['gene_id'].map(str)+"_"+longestDF['transcript_rank_in_gene'].astype(int).map(str)
    del(longestDF['transcript_rank_in_gene'])
    
    # Get transcript_per_gene for the consolidated set of transcripts
    longestDF['transcript_per_gene'] = longestDF.groupby('gene_id')['transcript_id'].transform('count')
    
    # Flag transcripts where the min start does not match the start of the longest representative
    longestDF['flag_not_min_start'] = np.where(longestDF['start']!=longestDF['min_transcript_start'],1,0)
    # Flag transcripts where the max end does not match the end of the longest representative
    longestDF['flag_not_max_end'] = np.where(longestDF['end']!=longestDF['max_transcript_end'],1,0)
    # Flag transcripts where a mono-exon gene was consolidated to pieces from multiple transcripts
    longestDF['flag_monoexon_collapse'] = np.where((longestDF['transcriptID_cat']==longestDF['junctionID_order'])&
             (longestDF['transcript_id'].str.contains("\|")),1,0)
    
    # Output key file for transcript_id 2 FSM_consolidation_transcript_id
    keyDF = FSM.split_transcript_id(df=longestDF.loc[:,['transcriptID_cat','FSM_consolidation_transcript_id']],
                                    col_name='transcriptID_cat')
    if len(keyDF) != len(xcrptGene):
        print("WARNING: Incorrect key file generated")
    del(xcrptGene)
    keyDF.to_csv(args.outDir+"/"+args.outPrefix+"_transcript_id_2_FSM_consolidation_transcript_id.csv",
                 index=False)
    
    # Print out counts for FSM consolidation
    print("{} transcripts and {} genes in FSM consolidation".format(longestDF['FSM_consolidation_transcript_id'].nunique(),longestDF['gene_id'].nunique()))
    
    # Output FSM consolidation transcript-level file
    longestDF.to_csv(args.outDir+"/"+args.outPrefix+"_FSM_consolidation.csv", index=False)

if __name__ == '__main__':
    # Parse command line arguments
    global args
    args = getOptions()
    main()

