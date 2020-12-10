#!/usr/bin/env python
# coding: utf-8

# In[327]:


import pandas as pd
import numpy as np
import argparse
import os
import matplotlib.pyplot as plt
import seaborn as sns
import math
from scipy.interpolate import make_interp_spline


# In[287]:


def getoptions():
    parser = argparse.ArgumentParser(description='calculate and visualize the distribution of feature number statistics from a standard file input')
    parser.add_argument("-ig", "--input_g", dest="input_g", required=True, help="Input file name")
    parser.add_argument("-ie", "--input_e", dest="input_e", required=True, help="Input file name")
    parser.add_argument("-if", "--input_f", dest="input_f", required=True, help="Input file name")
    parser.add_argument("-p", "--prefix", dest="prefix", required=True, help="a prefix for output file names")
    parser.add_argument("-G", "--columnG", dest="columnG", default='gene_id', help="specify the column name for unite feature")
    parser.add_argument("-T", "--columnT", dest="columnT", default='transcript_id', help="specify the column name for feature of which numbers will be counted")
    parser.add_argument("-F", "--columnF", dest="columnF", default='fragment_id', help="specify the column name for feature of which numbers will be counted")
    parser.add_argument("-E", "--columnE", dest="columnE", default='fusion_id', help="specify the column name for feature of which numbers will be counted")
    parser.add_argument("-o", "--output", dest="output", required=True, help="Output path")
    parser.add_argument("-s", "--separate", dest="separate", action='store_false', help="assign a feature to output a separate histgram figure")
    args = parser.parse_args()
    return(args)


# In[328]:


def remove_duprow(col, data):
    n = data[col].tolist()
    temp = []
    dup_index = []
    for i in range(0, len(n)):
        if n[i] not in temp:
            temp.append(n[i])
        else:
            dup_index.append(i)
    data = data.drop(index = dup_index)
    return data


# In[329]:


def calculate_number(unit, feature, data):
    nr = data.shape[0]
    temp = {}
    res = {}
    for i in range(0, nr):
        uni = data[unit][i].split("|")
        fea = data[feature][i].split("|")
        for j in uni:
            if j not in temp:
                temp[j] = fea
            else:
                for k in fea: 
                    if k not in temp[j]:
                        temp[j].append(k)
    for k in temp:
        res[k] = len(temp[k]) 
    return pd.Series(res)


# In[330]:


def calculate_num_from_pairfile(geneID, transID, annotDF):
    xcrptGene = annotDF[(~annotDF[transID].str.contains("|",regex=False))&
            (~annotDF[geneID].str.contains("|",regex=False))&
            (annotDF[transID]!="Unannotated")][[geneID, transID]].drop_duplicates()
    return xcrptGene.groupby(geneID).size()


# In[331]:


def generate_bins(data, nb = 10, merge_right = 99.5):
    bins = [min(data)]
    itv = (np.percentile(data, merge_right) - min(data))/(nb-1)
    for i in range(1, nb):
        bins.append(min(data) + i*itv)
    bins = bins + [max(data)]
    bins = [math.ceil(x) for x in bins]
    return bins


# In[374]:


def half_plot(dat, bins, colors, text=False):
    xs = []
    ys = []
    ylim = len(dat[dat.values <= bins[1]])* 1.1
    xlim = len(bins) * 1.1
    plt.xlim(1, xlim)
    plt.plot(xlim, ylim)
    plt.grid(False)
    plt.xticks([])
    for i in range(0, len(bins)-1):
        x = i+2
        y = len(dat[(dat.values < bins[i+1]) & (dat.values >= bins[i])])
        col = colors[i]
        plt.vlines(x, 0, y, color = col, linewidth=44)
        if text == True:
            plt.text(x, y, y, horizontalalignment='center',verticalalignment='bottom')
        xs.append(x)
        ys.append(y)
    return [xs, ys]


# In[394]:


def visualize_separate(list2, list3, out, prefix, kbins=10):
    titles = ['Number of exons regions per gene', 'Number of fragments per gene']
    names = ['_exonRegionsfig.pdf', '_fragmentfig.pdf']
    s = 0
    for i in [list2, list3]:
        bins = generate_bins(i[i>=2] , kbins, 99.9) # kbins (including a merged top 5% bin) + a bin for single transcript genes
        bins = [1] + bins
        fig = plt.figure(figsize=(10,6), facecolor='white',edgecolor='black')
        ax = fig.add_subplot(111)
        xs, ys = half_plot(i, bins, ['steelblue'] * (len(bins)-1), text=True)
        tick_pos = [1] + xs
        tick_pos = [i+0.5 for i in tick_pos]
        ax.set_xticks(tick_pos)
        ax.set_xticklabels([str(round(i, 2)) for i in bins] + [max(i)])
        ax.set_title(titles[s], fontsize=14)
        ax.grid(True)
        plt.savefig(out + '/' + prefix + names[s], dpi=600, format='pdf')
        s += 1


# In[379]:


def visualize_combined(list_main, list2, list3, out, prefix, text=True, kbins=10):
    dat = list_main
    dat2 = list2
    dat3 = list3
    xlabel = 'Number of transcripts / gene'
    ylabel = 'Number of genes'
    ylabel2 = 'Exon regions / gene'
    ylabel3 = 'Fragment / gene'
    titile = 'Distribution of number of transcripts per gene'

    bins = generate_bins(dat[dat>=2] , kbins, 99.9) # kbins (including a merged top 5% bin) + a bin for single transcript genes
    bins = [1] + bins

    font2 = {'family' : 'Times New Roman',
    'weight' : 'normal',
    'size'   : 15,
    }
    font3 = {'family' : 'Times New Roman',
    'weight' : 'normal',
    'size'   : 10,
    }

    ticklinewidth = 0.5
    colors = sns.color_palette("hls", kbins+1)
    height = len(dat[(dat.values < bins[1]) & (dat.values >= bins[0])]) * 1.1
    xlim = (kbins + 2) * 1.1

    sns.set_style("whitegrid")
    fig = plt.figure(figsize=(10,10), facecolor='white',edgecolor='black')
    grid = plt.GridSpec(6,1 , wspace=0, hspace=0.15)
    ys = []
    xs = []
    ###############################
    # plot section 1
    ax1 = fig.add_subplot(grid[0,0])
    ax1.patch.set_facecolor('white')
    ax1.axes.get_yaxis().set_visible(True)
    ax1.spines['bottom'].set_visible(False)
    ax1.set_title(titile, fontsize=14)

    xs, ys = half_plot(dat, bins, colors, text)
    
    break_y_up = np.percentile(ys, 90) * 0.7
    break_y_down = np.percentile(ys, 80) * 1.4
    ax1.set_ylim(break_y_up, height)
    tick_pos = [1] + xs
    tick_pos = [i+0.5 for i in tick_pos]

    ##################################
    # plot section 2
    ax2 = fig.add_subplot(grid[1:4,0])
    ax2.patch.set_facecolor('white')
    ax2.spines['top'].set_visible(False)
    ax2.spines['bottom'].set_visible(False)
    ax2.set_ylim(0, break_y_down)
    ax2.set_ylabel(ylabel,font2)
    ax2.plot(xs, ys, color = 'steelblue', lw = 1)
    #ys = []
    #xs = []
    half_plot(dat, bins, colors, text)

    ################################
    # extract genes based on bins
    feature2s = []
    feature3s = []
    for i in range(0, len(bins)-1):
        x = i+2
        genes = dat[(dat.values < bins[i+1]) & (dat.values >= bins[i])]
        feature2 = dat2[genes.index.tolist()].tolist()
        feature3 = dat3[genes.index.tolist()].tolist()
        feature2s.append(feature2)
        feature3s.append(feature3)

    #################
    # plot section 3
    ax3 = fig.add_subplot(grid[4,0])
    ax3.patch.set_facecolor('white')
    ax3.spines['top'].set_visible(False)
    ax3.spines['bottom'].set_visible(False)
    ax3.grid(False)
    ax3.set_ylim(0, max(dat2))
    ax3.set_ylabel(ylabel2, font3, labelpad = 11)
    box = ax3.boxplot(feature2s, positions = xs, notch = True, patch_artist=True)
    for patch, filers, whis, med, color in zip(box['boxes'], box['fliers'], box['whiskers'], box['medians'], colors):
        patch.set_color(color)
        filers.set_markeredgecolor(color)
        whis.set_color(color)
        med.set_color('white')

    ax3.set_xlim(1, xlim)
    ax3.set_xticklabels([])
    ax3.vlines(tick_pos, 0, max(dat2), color = 'grey', linewidth=ticklinewidth)

    y_arrow = max(dat2)
    for i in range(0, len(bins)-1): 
        ax3.arrow(xs[i], y_arrow, 0, -y_arrow/20, overhang=y_arrow/20, head_width=0.2, head_length=1, width = 1, shape="full",fc=colors[i], ec=colors[i],alpha=0.9)


    ##################################
    #plot section 4
    ax4 = fig.add_subplot(grid[5,0])
    ax4.patch.set_facecolor('white')
    ax4.spines['top'].set_visible(False)
    ax4.grid(False)
    ax4.set_ylim(0, max(dat3))
    ax4.set_ylabel(ylabel3, font3, labelpad = 6)
    box2 = ax4.boxplot(feature3s, positions = xs, notch = True, patch_artist=True)
    for patch, filers, whis, med, color in zip(box2['boxes'], box2['fliers'], box2['whiskers'], box2['medians'], colors):
        patch.set_color(color)
        filers.set_markeredgecolor(color)
        whis.set_color(color)
        med.set_color('white')

    y_arrow = max(dat3)
    ax4.set_xticks(tick_pos)
    ax4.set_xticklabels([str(round(i, 2)) for i in bins] + [max(dat)])
    ax4.set_xlim(1, xlim)
    ax4.set_xlabel(xlabel, font2, labelpad = 10)
    ax4.vlines(tick_pos, 0, max(dat3), color = 'grey', linewidth=ticklinewidth)
    for i in range(0, len(bins)-1): 
        ax4.arrow(xs[i], y_arrow, 0, -y_arrow/20, overhang=y_arrow/20, head_width=0.2, head_length=1, width = 1, shape="full",fc=colors[i], ec=colors[i], alpha=0.9)

    outfile = out + "/" + prefix + "_combinedfig.pdf" 
    plt.savefig(outfile, dpi=600, format='pdf')


# In[396]:


def main():
    args = getoptions()
    dat = pd.read_csv(args.input_g)
    num_list_main = calculate_num_from_pairfile(args.columnG, args.columnT, dat)
    del(dat)
    dat_e = pd.read_csv(args.input_e)
    #dat_e = remove_duprow(args.columnE, dat_e)
    num_list_e = calculate_number(args.columnG, args.columnE, dat_e)
    del(dat_e)
    dat_f = pd.read_csv(args.input_f)
    #dat_f = remove_duprow(args.columnF, dat_f)
    num_list_f = calculate_number(args.columnG, args.columnF, dat_f)
    del(dat_f)
    visualize_combined(num_list_main, num_list_e, num_list_f, out = os.path.normpath(args.output), prefix = args.prefix)
    if args.separate == True:
        visualize_separate(num_list_e, num_list_f, out = os.path.normpath(args.output), prefix = args.prefix)
    return


# In[ ]:


if __name__ == "__main__":
    main()



