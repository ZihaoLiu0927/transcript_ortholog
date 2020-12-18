#!/usr/bin/env python

#	DESCRIPTION: This program builds the database file required by GFFutils for GFF manipulation

## Modified by AVN
## Checks for gene and transcript features prior to database
##     generation to switch off gene and/or transcript inferring
##     and speed up database creation step


# Built-in packages
import argparse
import os

# Add-on packages
import gffutils

def getOptions():
    """ Function to pull in arguments """
    parser = argparse.ArgumentParser(description="Provide the path to the GFF file for creating the database")
    parser.add_argument("--gff", dest="gffInput", action='store', required=True, help="Input GFF file")
    
    args = parser.parse_args()
    return args

def main():
    gff_fn=args.gffInput
    db_fn=gff_fn + '.db'
    # check if db file exists and recreate if needed
    if os.path.isfile(db_fn):
        print("DB file exists. Deleting and recreating...")
        os.remove(db_fn)
    # Check if transcript or gene features are present in file
    disableTranscript = False
    disableGene = False
    gffFile = open(gff_fn, 'r')
    for line in gffFile:
        if len(line.split("\t")) >= 2:
            feature = line.split("\t")[2]
            if feature == "transcript":
                disableTranscript = True
            if feature == "gene":
                disableGene = True
            if disableTranscript and disableGene:
                break
        else:
            continue
    gffFile.close()
    gffutils.create_db(gff_fn, db_fn,merge_strategy='create_unique', disable_infer_transcripts=disableTranscript, disable_infer_genes=disableGene)
    #gffutils.create_db(gff_fn, db_fn,merge_strategy='create_unique')
    #gffutils.create_db(gff_fn, db_fn)

if __name__ == '__main__':
    # Parse command line arguments
    global args
    args = getOptions()
    main()

