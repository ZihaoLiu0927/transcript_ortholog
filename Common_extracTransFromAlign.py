import pandas as pd
import argparse
import os
import csv

def getoptions():
    parser = argparse.ArgumentParser(description='Extract counts information from lima outputs')
    parser.add_argument("-l", "--list", dest="list", required=True, help="file that contains a list of transcript names")
    parser.add_argument("-o", "--out", dest="out", required=True, help="Path to the directory where the output files will be in")
    parser.add_argument("-s", "--sam", dest="sam", required=True, help="A preffix for all the output files")
    args = parser.parse_args()
    return(args)


def parseSamHeader(sam):
    counts = 0
    with open(sam, "r") as f:
        for line in f:
            if line.startswith("@"):
                counts+=1
            else:
                return counts


def extractTranscript(list, sam, out):
    n_header = parseSamHeader(sam)
    transcript_list = pd.read_table(list, header = None, names = ["name"])
    samfile = pd.read_csv(sam, sep = "\n", header = None, names = ["name"])
    header_df = samfile.iloc[range(0, n_header),:]
    namecol = samfile.name.str.split("\t").str[0]
    newSAM = samfile[namecol.isin(transcript_list.name)]
    newSAM = header_df.append(newSAM)
    newSAM.to_csv(out, sep="\n", index = False, header = False, quoting=csv.QUOTE_NONE, quotechar = "", escapechar = '\n')
                    

def main():
    args = getoptions()
    extractTranscript(args.list, args.sam, args.out)
    return

if __name__ == "__main__":
    main()
