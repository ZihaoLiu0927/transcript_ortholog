#!/usr/bin/env python
# coding: utf-8


import pandas as pd
import numpy as np
import argparse
import matplotlib.pyplot as plt



def getoptions():
    parser = argparse.ArgumentParser(description='calculate and visualize the distribution of transcripts modified by different clustering methods')
    parser.add_argument('-a', "--c0mf", dest="c0mf", required=True, help="Input counts with label 0 by method 1")
    parser.add_argument('-b', "--c1mf", dest="c1mf", required=True, help="Input counts with label 1 by method 1")
    parser.add_argument('-c', "--c0mi", dest="c0mi", required=True, help="Input counts with label 0 by method 2")
    parser.add_argument('-d', "--c1mi", dest="c1mi", required=True, help="Input counts with label 1 by method 2")
    parser.add_argument('-o', "--outprefix", dest="outprefix", required=True, help="output file path and prefix")
    args = parser.parse_args()
    return(args)



def merge_fileCounts(a, b, c, d):
    df = a
    for i in [b,c,d]:
        df = pd.merge(df, i, how ="outer", on = "type")
    df.columns = ['type', 'flair_not_modified', "flair_modified", "isoseq3cluster_not_modified", "isoseq3cluster_modified"]
    df.index = df['type']
    df = df.fillna(0)
    df = df.drop('type', 1)
    return df



def generate_dfRatio(df):
    df_ratio = df.copy().fillna(0)
    for i in df_ratio.columns:
        df_ratio[i] = df_ratio[i]/df_ratio[i].sum()
    return df_ratio




def visualize_stackBar(df_ratio, outpdf = None, title = "Stack barplot of SQANTI classification ratio"):
    
    xs = [(i+1) for i in range(df_ratio.shape[1])]
    types = df_ratio.index
    height = 1
    ps = []
    colors = {'full-splice_match': "#FFC0CB", 
              'genic': "#87CEFA", 
              'incomplete-splice_match': "#7FFFAA", 
              'fusion': "#FF8C00", 
              'antisense': "#2F4F4F", 
              'novel_in_catalog': "#C0C0C0", 
              'novel_not_in_catalog': "#000000", 
              'intergenic': "#FF0000",
              'modified_both': "#FFC0CB", 
              'not_modified_both': "#87CEFA", 
              'modified_isoOnly': "#7FFFAA", 
              'modified_flairOnly': "#000000"}
    xtickslable = df_ratio.columns
    
    fig = plt.figure(figsize= (12,8))
    grid = plt.GridSpec(5,1 , wspace=0, hspace=0.15)
    ax = fig.add_subplot(grid[0:4, 0])
    
    xlims = (0,6.5 + df_ratio.shape[1] - 4)
    xticks = [i for i in range(1, df_ratio.shape[1]+1)]
    
    for i in types:
        p = ax.bar(xs, height, color = colors[i])
        ps.append(p)
        height = height - df_ratio.loc[i]
    ax.legend(ps, types, fontsize=10, loc=1) 
    ax.set_ylabel("Ratio", fontsize=14)
    ax.set_xlim(xlims)
    ax.set_title(title, fontsize=16)
    ax.set_xticks(xticks)
    ax.set_xticklabels(xtickslable, fontsize=14,rotation=35, horizontalalignment='right', verticalalignment='top')
    if not outpdf == None:
        fig.savefig(outpdf, pdi=600, format = 'pdf')
    plt.close(fig)


def main():
    args = getoptions()
    names = ["type", "counts"]
    a = pd.read_csv(args.c0mf, sep = "\t", names=names)
    b = pd.read_csv(args.c1mf, sep = "\t", names=names)
    c = pd.read_csv(args.c0mi, sep = "\t", names=names)
    d = pd.read_csv(args.c1mi, sep = "\t", names=names)
    
    df = merge_fileCounts(a, b, c, d)
    df_ratio = generate_dfRatio(df)
    
    outpdf = args.outprefix + "Stack_barplot_classRatio.pdf"
    visualize_stackBar(df_ratio, outpdf)

if __name__ == "__main__":
    main()




