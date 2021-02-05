#!/usr/bin/env python
# coding: utf-8


from FunctionsLibrary_ClusterStat import *


def getoptions():
    parser = argparse.ArgumentParser(description='count statistics for clustering methods')
    parser.add_argument('-f', "--flag", dest="flag", required=True, help="Standard input counts flag file")
    parser.add_argument('-c', "--cls", dest="cls", required=True, help="Input SQANTI QC classification file")
    parser.add_argument("-r", "--restrict", dest="restrict", default=None, help="a gene list file that will be used to subset the input data")
    parser.add_argument("-d", "--direction", dest="direction", default=None, help="select: more or less than the supporting read number condition")
    parser.add_argument('-o', "--outprefix", dest="outprefix", required=True, help="output file path and prefix")
    parser.add_argument('-g', "--groupby", dest="groupby", default = 'class', help="output file path and prefix")
    args = parser.parse_args()
    return(args)


def main():
    args = getoptions()
    
    dat = pd.read_csv(args.flag)
    cls = pd.read_csv(args.cls, sep = '\t', low_memory=False)
    cls = cls[["isoform", "structural_category"]]
    dat = dat.merge(cls, how = 'left', left_on= 'readID', right_on = 'isoform').drop('isoform', axis=1)
    readsCondition = "At least 2 "
    if not args.restrict == None:
        gene_filter = pd.read_csv(args.restrict, header = None, squeeze = True)
        if args.direction == "more":
            readsCondition = "At least 10 "
            dat = dat[dat['gene_id'].isin(gene_filter)].reset_index(drop=True)
        else:
            readsCondition = "less than 10 "
            dat = dat[~dat['gene_id'].isin(gene_filter)].reset_index(drop=True)
    clusterStat_table = make_geneCluster_numTable(dat)
    ###### output figures with different filtering criteria and different groupby methods
    
    groupby = args.groupby

    outGroupby_path = args.outprefix
    outGroupby_prefix = outGroupby_path + '_groupby_' + groupby + '_'
    i = 1


    outpdf = outGroupby_prefix + str(i) + '.pdf'
    count, ratio, ratio_removed = merge_Counts(dat, criteria = 'OneCluster', check_FSM = None, groupby = groupby)
    title = readsCondition + "reads, one cluster per gene in either method"
    if not count.empty:
        visualize_stackBar(ratio, count, outpdf = outpdf, title = title)
        count.to_csv(outpdf.replace(".pdf", '.csv'))
        ratio_removed.round(decimals = 2).to_csv(outpdf.replace(".pdf", '_removedRatio.csv'))
    else:
        print("condition without output: " + title)
    i+=1


    outpdf = outGroupby_prefix + str(i) + '.pdf'
    count, ratio, ratio_removed = merge_Counts(dat, criteria = 'OneCluster', check_FSM = 'AtLeastOneEither', groupby = groupby)
    title = readsCondition + "reads, one cluster per gene in either method, at least 1 FSM"
    if not count.empty:
        visualize_stackBar(ratio, count, outpdf = outpdf,title = title)
        count.to_csv(outpdf.replace(".pdf", '.csv'))
        ratio_removed.round(decimals = 2).to_csv(outpdf.replace(".pdf", '_removedRatio.csv'))
    else:
        print("condition without output: " + title)
    i+=1


    outpdf = outGroupby_prefix + str(i) + '.pdf'
    count, ratio, ratio_removed = merge_Counts(dat, criteria = 'MoreThanOne', check_FSM = None, groupby = groupby)
    title = readsCondition + "reads, at least 2 clusters per gene in either method"    
    if not count.empty:
        visualize_stackBar(ratio, count, outpdf = outpdf, title = title)
        count.to_csv(outpdf.replace(".pdf", '.csv'))
        ratio_removed.round(decimals = 2).to_csv(outpdf.replace(".pdf", '_removedRatio.csv'))
    else:
        print("condition without output: " + title)
    i+=1


    outpdf = outGroupby_prefix + str(i) + '.pdf'
    count, ratio, ratio_removed = merge_Counts(dat, criteria = 'MoreThanOne', check_FSM = 'AtLeastOneEither', groupby = groupby)
    title = readsCondition + "reads, at least 2 clusters per gene in either method, at least 1 FSM"    
    if not count.empty:
        visualize_stackBar(ratio, count, outpdf = outpdf, title = title)
        count.to_csv(outpdf.replace(".pdf", '.csv'))
        ratio_removed.round(decimals = 2).to_csv(outpdf.replace(".pdf", '_removedRatio.csv'))
    else:
        print("condition without output: " + title)
    i+=1


    outpdf = outGroupby_prefix + str(i) + '.pdf'
    count, ratio, ratio_removed = merge_Counts(dat, criteria = 'MoreThanOne', check_FSM = 'AtLeastOneBoth', groupby = groupby)
    title = readsCondition + "reads, at least 2 clusters per gene in either method, at least 1 FSM in both cluster"
    if not count.empty:
        visualize_stackBar(ratio, count, outpdf = outpdf, title = title)
        count.to_csv(outpdf.replace(".pdf", '.csv'))
        ratio_removed.round(decimals = 2).to_csv(outpdf.replace(".pdf", '_removedRatio.csv'))
    else:
        print("condition without output: " + title)
    i+=1


    clusterStat_table = make_geneCluster_numTable(dat)


    outpdf = outGroupby_prefix + str(i) + '.pdf'
    count, ratio, ratio_removed = merge_Counts(dat, criteria = 'EqualOneCluster', clusterStat = clusterStat_table, check_FSM = None, groupby = groupby)
    title = readsCondition + "reads, at least 1 equal cluster per gene in both methods"
    if not count.empty:
        visualize_stackBar(ratio, count, outpdf = outpdf, title = "At least 2 reads, 1 equal cluster per gene in both methods")
        count.to_csv(outpdf.replace(".pdf", '.csv'))
        ratio_removed.round(decimals = 2).to_csv(outpdf.replace(".pdf", '_removedRatio.csv'))
    else:
        print("condition without output: " + title)
    i+=1


    outpdf = outGroupby_prefix + str(i) + '.pdf'
    count, ratio, ratio_removed = merge_Counts(dat, criteria = 'EqualOneCluster', clusterStat = clusterStat_table, check_FSM = 'AtLeastOneBoth', groupby = groupby)
    title = readsCondition + "reads, 1 equal cluster per gene in both methods, at least 1 FSM in each cluster"
    if not count.empty:
        visualize_stackBar(ratio, count, outpdf = outpdf, title = title)
        count.to_csv(outpdf.replace(".pdf", '.csv'))
        ratio_removed.round(decimals = 2).to_csv(outpdf.replace(".pdf", '_removedRatio.csv'))
    else:
        print("condition without output: " + title)
    i+=1


    outpdf = outGroupby_prefix + str(i) + '.pdf'
    count, ratio, ratio_removed = merge_Counts(dat, criteria = 'EqualOneCluster', clusterStat = clusterStat_table, check_FSM = 'NoOneNeither', groupby = groupby)
    title = readsCondition + "reads, 1 equal cluster per gene in both methods, without FSM in each cluster"
    if not count.empty:
        visualize_stackBar(ratio, count, outpdf = outpdf, title = title)
        count.to_csv(outpdf.replace(".pdf", '.csv'))
        ratio_removed.round(decimals = 2).to_csv(outpdf.replace(".pdf", '_removedRatio.csv'))
    else:
        print("condition without output: " + title)
    i+=1


    outpdf = outGroupby_prefix + str(i) + '.pdf'
    count, ratio, ratio_removed = merge_Counts(dat, criteria = 'EqualMoreThanOne', clusterStat = clusterStat_table, check_FSM = None, groupby = groupby)
    title = readsCondition + "reads, at least 2 equal clusters per gene in both methods"
    if not count.empty:
        visualize_stackBar(ratio, count, outpdf = outpdf, title = title)
        count.to_csv(outpdf.replace(".pdf", '.csv'))
        ratio_removed.round(decimals = 2).to_csv(outpdf.replace(".pdf", '_removedRatio.csv'))
    else:
        print("condition without output: " + title)
    i+=1


    outpdf = outGroupby_prefix + str(i) + '.pdf'
    count, ratio, ratio_removed = merge_Counts(dat, criteria = 'EqualMoreThanOne', clusterStat = clusterStat_table, check_FSM = 'AtLeastOneBoth', groupby = groupby)
    title = readsCondition + "reads, at least 2 equal clusters per gene in both methods, at least 1 FSM in each cluster"
    if not count.empty:
        visualize_stackBar(ratio, count, outpdf = outpdf, title = title)
        count.to_csv(outpdf.replace(".pdf", '.csv'))
        ratio_removed.round(decimals = 2).to_csv(outpdf.replace(".pdf", '_removedRatio.csv'))
    else:
        print("condition without output: " + title)
    i+=1


    outpdf = outGroupby_prefix + str(i) + '.pdf'
    count, ratio, ratio_removed = merge_Counts(dat, criteria = 'EqualMoreThanOne', clusterStat = clusterStat_table, check_FSM = 'NoOneNeither', groupby = groupby)
    title = readsCondition + "reads, at least 2 equal clusters per gene in both methods, without FSM in each cluster"
    if not count.empty:
        visualize_stackBar(ratio, count, outpdf = outpdf, title = title)
        count.to_csv(outpdf.replace(".pdf", '.csv'))
        ratio_removed.round(decimals = 2).to_csv(outpdf.replace(".pdf", '_removedRatio.csv'))
    else:
        print("condition without output: " + title)


if __name__ == '__main__':
    main()
