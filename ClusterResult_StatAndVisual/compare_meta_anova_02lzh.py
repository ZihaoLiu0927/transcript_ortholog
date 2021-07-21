#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun  3 13:16:50 2021

@author: zihaoliu
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib.patches as patches
import argparse

def getOptions():
    parser = argparse.ArgumentParser(description="get parameters from Galaxy")
    parser.add_argument('-m', '--model', dest='model', required=True, help="The ANOVA model")
    parser.add_argument('-e', '--efm', dest='efm', required=True, help="meta-analysis effect size")
    parser.add_argument('-r', '--rank', dest='rank', default=None, help="The rank method, can be logPH, rank1000 or rank100")
    parser.add_argument('-s', '--sig', dest='sig', default=0.05, type = float, help="significant shreshod")
    args=parser.parse_args()
    return args



def plot_pie(ax, counts, sig_group, outer_colors):
    size=0.45
    patches, texts = ax.pie(x = counts, labels =  
                None, labeldistance = 0.5, shadow = False, colors = outer_colors,
                        radius = 1.22, wedgeprops=dict(width=size, edgecolor='w', linewidth=0.2))
        
    #patches2, texts2 = ax.pie(x = df_counts_opposite['all_met'], labels =  
    #            None, labeldistance = 0.2, 
    #            colors = inner_colors, startangle = 180, shadow = False,
    #            radius=1-size, wedgeprops=dict(width=size, edgecolor='w'))

    legend1 = plt.legend(patches, sig_group, loc='center', prop={'size': 5},  framealpha = 0.3, frameon=False)
    #legend2 = plt.legend(patches2, df_counts_opposite['all_met'], loc='upper left', prop={'size': 6}, framealpha = 0.5)
    ax.add_artist(legend1)
    #ax.add_artist(legend2)
    


def plot_scatter(ax, x, y, colors, labels, xlabel, ylabel, pos_vline = None, pos_hline = None, rec=False, legend = False, xlim = None, ylim = None, pointSize = [1.2, 1.2], pshape = ['1', '1'], color_labels = 'colors', zorder = 3):
    if xlim == None:
        xlim = (min(x), max(x))
    if ylim == None:
        ylim = (min(y), max(y))
    
    for i in data[color_labels].unique():
        idx = data[color_labels] == i
        ax.plot(x[idx], y[idx], pshape[i], color = colors[i], markersize = pointSize[i], label=labels[i])
    if pos_vline:
        ax.vlines(pos_vline, ylim[0], ylim[1], 'r', linestyle = 'dashed', linewidths = 0.8, zorder=zorder)
    if pos_hline:
        ax.hlines(pos_hline, xlim[0], xlim[1], 'r', linestyle = 'dashed', linewidths = 0.8, zorder=zorder)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    

    xticks = list(map(round, [xlim[0], xlim[0]/2, 0, xlim[1]/2 , xlim[1]]))
    #yticks = list(map(round, [ylim[0], ylim[0]/2 ,0, ylim[1]/2 , ylim[1]]))
    rotation =0
    markerscale = 2
    
    ax.set_xticks(xticks)
    #ax.set_yticks(yticks)
    ax.set_xticklabels(xticks, rotation=rotation)
    #ax.set_yticklabels(yticks)

    if rec:
        rect = patches.Rectangle((0, 0), 0.05, 0.05, linewidth=1.5, edgecolor=None, facecolor='red', alpha=0.5)
        ax.add_patch(rect)
    if legend:
        ax.legend(loc='upper left', markerscale = markerscale, framealpha = 0.3, prop={'size': 7})
    return ax


def myplot(out):
    fig = plt.figure(figsize=(8,8), facecolor='white',edgecolor='black')
    grid = plt.GridSpec(7, 7, wspace=0, hspace=0)
    ax_scatter = fig.add_subplot(grid[2:7,0:5])
    ax_x_1 = fig.add_subplot(grid[1, 0:5])
    ax_x_2 = fig.add_subplot(grid[0, 0:5])
    ax_y_1 = fig.add_subplot(grid[2: 7,5])
    ax_y_2 = fig.add_subplot(grid[2: 7,6])
    ax_stat_all = fig.add_subplot(grid[0:2, 5:7])
    ax_x_1.get_xaxis().set_visible(False)
    ax_x_1.get_yaxis().set_visible(False)
    ax_y_1.get_xaxis().set_visible(False)
    ax_y_1.get_yaxis().set_visible(False)
    ax_x_2.get_xaxis().set_visible(False)
    ax_y_2.get_yaxis().set_visible(False)
    
    #plot_pie(ax_stat_all, sig_group, df_counts_opposite, sigGroup_color, colors)
    
    ax_scatter = plot_scatter(ax_scatter, x, y, group_colors, group_legends, model_x, model_y, pos_vline = 0.0001,
                      rec=False, legend=True, pointSize = [2] * 5, pshape = ['o','o', 'o', 'o'])
    
    xlabel_sub = "p sig in meta" 
    ylabel_sub = "p sig in " + model
    
    #bins_y = np.arange(up_show, bot_show, (abs(up_show) + abs(bot_show))/25)
    #bins_x = np.arange(left_show, right_show, (abs(left_show) + abs(right_show))/25)  
    
    #sns.histplot(ax = ax_x_1, x = x, kde=True, stat = 'count', edgecolor = '#E0E0E0', bins = bins_x)
    sns.histplot(ax = ax_x_1, x = x, kde=True, stat = 'count', edgecolor = '#E0E0E0')
    plot_scatter(ax_x_2, px, py, group_colors, group_legends, xlabel_sub, ylabel_sub, pos_vline = 0.05, xlim = (0, 1), ylim = (0, 0.05), pointSize=[2] * 5, pshape = ['o'] * 5)
    #sns.histplot(ax = ax_y_1, y = y, kde=True, stat = 'count', edgecolor = '#E0E0E0', bins = bins_y)
    sns.histplot(ax = ax_y_1, y = y, kde=True, stat = 'count', edgecolor = '#E0E0E0')
    plot_scatter(ax_y_2, px, py, group_colors, group_legends, xlabel_sub, ylabel_sub, pos_hline = 0.05, xlim = (0, 0.05), ylim = (0, 1), pointSize=[2] * 5, pshape = ['o'] * 5)
    plot_pie(ax_stat_all, counts, group_legends, group_colors[0:3] + ["black"])
    fig.set_rasterized(True)
    fig.savefig(out, dpi=300, format = 'pdf')


def main():
    args = getOptions()
    model_names = {'CtISB' : 'Contrast in a set and batch covariate',
                   'CtIS' : 'Contrast in a set',
                   'CtISP' : 'Contrast in a set and pairing covariate',
                   'pr' : 'anova form paired t-test',
                   'unpr' : 'regular anova',
                   'CISB' : 'control in a set',
                   'CISP' : 'control in a set and batch covariate',
                   'CtAS' : "contrast in all set and batch covariate",
                   'CtASB' : 'contrast in all set and batch covariate',
                   'CtASP' : 'contrast in all set and batch covariate',
    }

    mutinfo = pd.read_table("~/mclab/SHARE/McIntyre_Lab/CID/sweet16/design_files/batch_model_tests/all_mutantsList.tsv", header = None)
    muts = mutinfo[0].tolist()
    model = args.model
    efm = args.efm
    alpha = args.sig
    rank = args.rank

    for mut in muts:
        print(mut)
        file_path = "~/mclab/SHARE/McIntyre_Lab/CID/sweet16/models_1442/all_model_tests/model_output/merged_result/drop_miss_std_" + mut + "_" + model +"_" + efm + "_" + rank +  "_merged_anova_meta.tsv"
        data = pd.read_csv(file_path, sep = '\t')
        
        xlab = 'meta_effect'
        ylab = 'StdMD'
        #ylab = 'Estimate'
        
        model_x = 'meta analysis effect size'
        model_y = model + " lsmean difference diviced by standard deviation"
    
        
        xplab = 'meta_p_value'
        yplab = 'pv_' + model + "_" + mut
        px = data[xplab]
        py = data[yplab]
        ################
        data_xy_both = data[(px < alpha) & (py < alpha)]
        data_x_only = data[(px < alpha) & (py >= alpha)]
        data_y_only = data[(px >= alpha) & (py < alpha)]
        data_xy_neither = data[(px > alpha) & (py > alpha)]
        counts = [len(x) for x in [data_xy_both, data_x_only, data_y_only, data_xy_neither]]
        
    
        group_colors = ['#FF0000', '#3399FF','#FF9933', '#808080']
        
        lgx = "sig_in_meta_only"
        lgy = "sig_in_" + model + "_only"
        group_legends = ['sig_both', lgx, lgy, 'sig_neither']
        group_legends = [ (x + ": " +str(y)) for x, y in zip(group_legends, counts)]
        #######
        
        data["colors"] = 0
        data.loc[data_x_only.index, 'colors'] = 1
        data.loc[data_y_only.index, 'colors'] = 2
        data.loc[data_xy_neither.index, 'colors'] = 3
                  
        x = data[xlab]
        y = data[ylab]

        outfile = "/nfshome/zihaoliu/mclab/SHARE/McIntyre_Lab/CID/sweet16/models_1442/all_model_tests/model_output/compare_models/ANOVA_meta/drop_miss_std_" + mut + "_" + model +  "_" + efm + "_" + rank +  "_scatterplot.pdf" 
        myplot(out = outfile)


if __name__ == "__main__":
    main()
