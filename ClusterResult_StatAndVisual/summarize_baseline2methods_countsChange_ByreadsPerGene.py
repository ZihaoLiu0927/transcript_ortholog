#!/usr/bin/env python
# coding: utf-8


from FunctionsLibrary_ClusterStat import *



def getoptions():
    parser = argparse.ArgumentParser(description='count statistics for clustering methods')
    parser.add_argument('-f', "--flag", dest="flag", required=True, help="Standard input counts flag file")
    parser.add_argument('-c', "--cls", dest="cls", required=True, help="Input SQANTI QC classification file")
    parser.add_argument('-o', "--outprefix", dest="outprefix", required=True, help="output file path and prefix")
    parser.add_argument('-a', "--rangeA", dest="rangeA", required=True, help="supporting reads range from")
    parser.add_argument('-b', "--rangeB", dest="rangeB", required=True, help="supporting reads range to")
    args = parser.parse_args()
    return(args)


def main():
    args = getoptions()
    dat = pd.read_csv(args.flag)
    cls = pd.read_csv(args.cls, sep = '\t', low_memory=False)
    cls = cls[["isoform", "structural_category"]]
    dat = dat.merge(cls, how = 'left', left_on= 'readID', right_on = 'isoform').drop('isoform', axis=1)
    df_sum = make_df_clusterCompare_BysupReads(int(args.rangeA), int(args.rangeB), dat)
    outfile = args.outprefix + '_geneWithN_sptReads_ClusterPerGene_counts.txt'
    df_sum.to_csv(outfile)
    return

if __name__ == '__main__':
    main()




