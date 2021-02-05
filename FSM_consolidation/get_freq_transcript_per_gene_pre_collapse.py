#!/usr/bin/env python

import argparse
import pandas as pd

def getOptions():
    # Parse command line arguments
    parser = argparse.ArgumentParser(description="Get frequency of transcripts per gene prior to collapse")

    # Input data
    parser.add_argument("-i", "--input", dest="inXcrptGene", required=True, help="CSV file with unique pairs of transcript_id to gene_id")

    # Output data
    parser.add_argument("-o", "--output", dest="outFile", required=True, help="Output CSV file with the frequencies of transcripts per gene")

    args = parser.parse_args()
    return args

def main():
    # Get input file of event_id to transcript_id to gene_id
    xcrptGene = pd.read_csv(args.inXcrptGene, low_memory=False)

    # Output total number of transcripts and genes pre-collapse
    print("\n{} genes pre-collapse\n{} transcripts pre-collapse\n".format(
            xcrptGene['gene_id'].nunique(),xcrptGene['transcript_id'].nunique()))

    # Add transcripts per gene variable
    xcrptGene['num_transcript_per_gene'] = xcrptGene.groupby('gene_id')['transcript_id'].transform('count')
    
    # Get frequency of transcripts per gene prior to collapse
    xcrptGene.groupby('num_transcript_per_gene')['gene_id'].nunique().to_csv(
            args.outFile, index_label = 'transcript_per_gene', header=['freq'])
    
if __name__ == '__main__':
    # Parse command line arguments
    global args
    args = getOptions()
    main()

