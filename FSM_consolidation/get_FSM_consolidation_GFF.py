#!/usr/bin/env python

import argparse
import gffutils
import pandas as pd

# Import custon functions for FSM consolidation
import FSM_consolidation_functions as FSM

def getOptions():
    # Parse command line arguments
    parser = argparse.ArgumentParser(description="Generate FSM consolidated GFF files for single-transcript and multi-transcript genes (output files are unsorted)")

    # Input data
    parser.add_argument("-i", "--input-FSM-consolidation", dest="inFSM", required=True, help="CSV file from FSM consolidation")
    parser.add_argument("-g", "--gff-db", dest="inDB", required=True, help="Converted GFF database file (*.gff.db) created in event analysis annotation generation")

    # Output data
    parser.add_argument("-d", "--output-directory", dest="outDir", required=True, help="Output directory")
    parser.add_argument("-p", "--output-prefix", dest="outPrefix", required=True, help="Prefix for output files")

    args = parser.parse_args()
    return args

def main():
    # Get FSM consolidation input file
    longestDF = pd.read_csv(args.inFSM, low_memory=False)

    # Get coverted GFF database from event analysis annotation generation
    db = gffutils.FeatureDB(args.inDB)
    
    # Open output GTF files
    singleGTF = open(args.outDir+"/"+args.outPrefix+"_FSM_consolidation_single_transcript.gff", 'w')
    multiGTF = open(args.outDir+"/"+args.outPrefix+"_FSM_consolidation_multi_transcript.gff", 'w')
    
    # Group FSM consolidation transcripts by gene_id
    genes = longestDF.groupby('gene_id')
    
    # Loop over all genes in FSM consolidation transcriptome
    for name, group in genes:
        # Check if gnee is single-transcript or multi-transcript gene
        #     and set proper output file
        if group['transcript_per_gene'].unique()[0] == 1:
            outFile = singleGTF
        else:
            outFile = multiGTF
        # Output gene feature
        feature = db[name]
        feature.source = "FSM_consolidation"
        outFile.write(str(feature)+"\n")
        # Set dataframe for exons in gene
        exonDF = pd.DataFrame(columns=['chr','source','feature','start','end',
                                         'score','strand','phase','attribute_Name',
                                         'attribute_Parent','attribute_parent_type'])
        # Loop over all FSM consolidation transcripts within each gene,
        #     modify where needed, and output
        for index, row in group.iterrows():
            # Get transcript feature (will take the first value in piped list of monoexon transcript_id values)
            transcript = row['transcript_id'].split("|")[0]
            feature = db[transcript]
            if row['flag_not_min_start']==1:
                if row['flag_not_max_end']==1:
                    feature = FSM.modify_gff_attributes(feature,row,diff_start=row['min_transcript_start'],
                                                    diff_end=row['max_transcript_end'])
                else:
                    feature = FSM.modify_gff_attributes(feature,row,diff_start=row['min_transcript_start'])
            else:
                if row['flag_not_max_end']==1:
                    feature = FSM.modify_gff_attributes(feature,row,diff_end=row['max_transcript_end'])
                else:
                    feature = FSM.modify_gff_attributes(feature,row)
            if row['flag_monoexon_collapse']==1:
                feature = FSM.modify_gff_attributes(feature,row,diff_start=row['min_transcript_start'],
                                                    diff_end=row['max_transcript_end'])
            outFile.write(str(feature)+"\n")
            
            # Loop over exons of transcript, modify gff attributes/coordinates, and output to file
            count = 0
            total = len(list(db.children(feature, featuretype='exon')))
            for exon in db.children(feature, featuretype='exon', order_by='start'):
                count = count + 1
                # Mono-exon collapse transcripts will need the exon modified on both ends
                if row['flag_monoexon_collapse']==1 and total == 1:
                    exon = FSM.modify_gff_attributes(exon,row,diff_start=row['min_transcript_start'],
                                                    diff_end=row['max_transcript_end'])
                # Modify first exon of transcripts with modified 5' end
                elif row['flag_not_min_start']==1 and count==1:
                    exon = FSM.modify_gff_attributes(exon, row, diff_start=row['min_transcript_start'])
                # Modify last exon of transcripts with modified 3' end
                elif row['flag_not_max_end']==1 and count==total:
                    exon = FSM.modify_gff_attributes(exon, row, diff_end=row['max_transcript_end'])
                else:
                    exon = FSM.modify_gff_attributes(exon, row)
                # Check if exon is already present in another transcript
                # If yes then add new transcript to Parent list
                # If no then add exon to exon list
                exonDF = FSM.add_exon(exon, exonDF)
                
            # Loop over introns of transcript modify Parent, and output to file
            # (none are modified since internal splice sites are not affected in consolidation)
            for intron in db.children(feature, featuretype='intron', order_by='start'):
                intron = FSM.modify_gff_attributes(intron, row)
                outFile.write(str(intron)+"\n")
                
        # Print all exons in gene to output file
        for index, row in exonDF.iterrows():
            outFile.write("\t".join(row[['chr','source','feature','start','end',
                        'score','strand','phase']].astype(str))+"\t"+";".join(
                        row[['attribute_Name','attribute_Parent','attribute_parent_type']].astype(str))+"\n")
    singleGTF.close()
    multiGTF.close()
    
if __name__ == '__main__':
    # Parse command line arguments
    global args
    args = getOptions()
    main()

