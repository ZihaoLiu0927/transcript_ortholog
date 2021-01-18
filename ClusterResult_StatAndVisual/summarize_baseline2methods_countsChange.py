#!/usr/bin/env python
# coding: utf-8

from FunctionsLibrary_ClusterStat import *


def getoptions():
    parser = argparse.ArgumentParser(description='count statistics for clustering methods')
    parser.add_argument('-f', "--flag", dest="flag", required=True, help="Standard input counts flag file")
    parser.add_argument('-c', "--cls", dest="cls", required=True, help="Input SQANTI QC classification file")
    parser.add_argument('-o', "--outprefix", dest="outprefix", required=True, help="output file path and prefix")
    args = parser.parse_args()
    return(args)


def main():
    args = getoptions()
    
    dat = pd.read_csv(args.flag)
    cls = pd.read_csv(args.cls, sep = '\t', low_memory=False)
    cls = cls[["isoform", "structural_category"]]
    dat = dat.merge(cls, how = 'left', left_on= 'readID', right_on = 'isoform').drop('isoform', axis=1)
    dat_sr = dat[dat['flag_single_read_gene_baseline'] == 1].reset_index()
    
    single_read_genes_table = single_read_genes_Stat(dat_sr)
    single_read_gene_cls_table = single_read_gene_cls_Stat(dat_sr)
    
    clusterStat_table = make_geneCluster_numTable(dat)
    clusterStat_out = clusterStat_table[[i for i in clusterStat_table.columns if 'flag' in i]].sum().to_frame()
    clusterStat_out.columns = ['counts']
    
    dat_equal = dat[(dat['gene_id'].isin(clusterStat_table[clusterStat_table['flag_equalNum'] == 1].index))]
    equal_list_iso = dat_equal[(dat_equal['flag_multigene_isoseq3Cluster'] == 0) & 
                (dat_equal['flag_in_isoseq3Cluster'] == 1)].groupby('gene_id')['iso_cluster_id'].unique().apply(
                lambda x : '|'.join(x)).reset_index()
    equal_list_flair = dat_equal[(dat_equal['flag_multigene_FLAIR'] == 0) & 
                (dat_equal['flag_in_FLAIR'] == 1)].groupby('gene_id')['flair_transcript_id'].unique().apply(
                lambda x : '|'.join(x)).reset_index()
    equalGeneList = equal_list_iso.merge(equal_list_flair, on = 'gene_id')
    
    outfile_sr = args.outprefix  + '_Single_read_gene_counts.txt'
    outfile_clusterCount = args.outprefix  + '_ClusterPerGene_counts.txt'
    outfile_geneList = args.outprefix  + '_EqualCluster_geneList.txt'
    #########
    #### output counts tables
    equalGeneList.to_csv(outfile_geneList, index = False)
    with open(outfile_sr, 'w') as f:
        f.write("Single read genes counts: \n")
        f.write("{}\n".format(single_read_genes_table.to_string()))
        f.write("\n\nSQANTI classification of the single read genes removed: \n")
        f.write("{}\n".format(single_read_gene_cls_table.to_string()))
    with open(outfile_clusterCount, 'w') as f:
        f.write("Cluster per gene comparison counts: \n")
        f.write("{}\n".format(clusterStat_out.to_string()))
    ##########


if __name__ == '__main__':
    main()





