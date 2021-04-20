#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 29 16:40:28 2021

@author: zihao liu
"""

import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np
import seaborn as sns



types = 3
md_col_first = 'diff_of_PD1074-RB2011'
md_col_second = 'diff_of_PD1074-RB2011.1'
pvalue_col_first = 'prob_greater_than_t_for_diff_PD1074-RB2011'
pvalue_col_second = 'prob_greater_than_t_for_diff_PD1074-RB2011.1'
tvalue_direct_first = 't-Value_for_Diff_PD1074-RB2011'
tvalue_direct_second = 't-Value_for_Diff_PD1074-RB2011.1'
xlabel = 'anova raw scale mean_difference'
ylabel = 'anova rank bin 1000 mean_difference'
file = "/Users/zach/Desktop/anovaRowScale_anovaRankBin1000_cmp_1.pdf"

data = pd.read_table("/Users/zach/Desktop/anovaRowScale_anovaRankBin1000_merge_set1_HILIC_pos_RB2011.tsv")
data['opposite_direction'] = np.where(data[tvalue_direct_first] * data[tvalue_direct_second] > 0, 0, 1)
data['color_labels'] = data['opposite_direction'] 



df_counts_opposite = pd.DataFrame({'all_met' : [0, 0]})
df_counts_opposite.index = ['same_direction', 'opposite_direction']
dic_direct = {0: "same_direction", 1: "opposite_direction"}
counts = data['opposite_direction'].value_counts()
counts.index = list(map(lambda x : dic_direct[x], counts.index))
for i in counts.index:
    df_counts_opposite.loc[i] = counts[i]
df_counts_opposite = pd.DataFrame(df_counts_opposite)
df_counts_opposite.columns = ['all_met']



data_cutoff_both = data[(data[pvalue_col_first] < 0.05) & (data[pvalue_col_second] < 0.05)]
data_cutoff_xaxis = data[(data[pvalue_col_first] < 0.05) & (data[pvalue_col_second] >= 0.05)]
data_cutoff_yaxis = data[(data[pvalue_col_first] >= 0.05) & (data[pvalue_col_second] < 0.05)]
data_cutoff_neither = data[(data[pvalue_col_first] > 0.05) & (data[pvalue_col_second] > 0.05)]
sig_group = pd.DataFrame({"sig_both" : len(data_cutoff_both), 
                          "sig_xaxis_only" : len(data_cutoff_xaxis), 
                          "sig_yaxis_only" : len(data_cutoff_yaxis), 
                          "sig_neither" : len(data_cutoff_neither)}, 
                         index=["counts"]).T

group_colors = ['#3399FF', '#FF9933', '#FF0000', '#00FF00', '#9370DB']
group_labels = ['same_direction', 'opposite_direction', 'sig_both', 'sig_xaxis_only', 'sig_yaxis_only']
group_shape = ['1','1','o', 'o', 'o']

if types != 3:
    x = data[pvalue_col_first]
    y = data[pvalue_col_second]
    x_same = data[data['opposite_direction'] == 0][pvalue_col_first]
    y_same = data[data['opposite_direction'] == 0][pvalue_col_second]
    x_oppo = data[data['opposite_direction'] == 1][pvalue_col_first]
    y_oppo = data[data['opposite_direction'] == 1][pvalue_col_second]
    
    colors = group_colors[0:2]
    labels = group_labels[0:2]
    pshape = group_shape[0:2]
    
else:
    x = data[md_col_first]
    y = data[md_col_second]
    x_same = data[data['opposite_direction'] == 0][md_col_first]
    y_same = data[data['opposite_direction'] == 0][md_col_second]
    x_oppo = data[data['opposite_direction'] == 1][md_col_first]
    y_oppo = data[data['opposite_direction'] == 1][md_col_second]
    
    
    data.loc[data_cutoff_both.index, 'color_labels'] = 2
    data.loc[data_cutoff_xaxis.index, 'color_labels'] = 3
    data.loc[data_cutoff_yaxis.index, 'color_labels'] = 4
    
    colors = group_colors
    labels = group_labels
    pshape = group_shape




def plot_scatter(ax, x, y, data, colors, labels, xlabel, ylabel, pos_vline, pos_hline, rec=False, legend = False, xlim = None, ylim = None, pointSize = [1.2, 1.2], pshape = ['1', '1'], types = 1):
    if xlim == None:
        xlim = (min(x), max(x))
    if ylim == None:
        ylim = (min(y), max(y))
    
    for i in data['color_labels'].unique():
        idx = data['color_labels'] == i
        ax.plot(x[idx], y[idx], pshape[i], color = colors[i], markersize = pointSize[i], label=labels[i])
        
    ax.vlines(pos_vline, ylim[0], ylim[1], 'black', linestyle = 'dashed', linewidths = 0.7, zorder=3)
    ax.hlines(pos_hline, xlim[0], xlim[1], 'black', linestyle = 'dashed', linewidths = 0.7, zorder=3)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    
    if types != 3:
        xticks = [0, xlim[1]/2 , xlim[1]]
        yticks = [0, ylim[1]/2 , ylim[1]]
        rotation = 0
        markerscale = 4
    else:
        xticks = list(map(round, [xlim[0], xlim[0]/2, 0, xlim[1]/2 , xlim[1]]))
        yticks = list(map(round, [ylim[0], ylim[0]/2 ,0, ylim[1]/2 , ylim[1]]))
        rotation =45
        markerscale = 2
    
    ax.set_xticks(xticks)
    ax.set_yticks(yticks)
    ax.set_xticklabels(xticks, rotation=rotation)
    ax.set_yticklabels(yticks)

    if rec:
        rect = patches.Rectangle((0, 0), 0.05, 0.05, linewidth=1.5, edgecolor=None, facecolor='red', alpha=0.5)
        ax.add_patch(rect)
    if legend:
        ax.legend(loc='upper left', markerscale = markerscale, framealpha = 0.3, prop={'size': 7})
    return ax


def plot_pie(ax, sig_group, df_counts_opposite, outer_colors, inner_colors):
    size=0.4
    patches, texts = ax.pie(x = sig_group['counts'], labels =  
                None, labeldistance = 0.5, shadow = False, colors = outer_colors,
                radius=1, wedgeprops=dict(width=size, edgecolor='w'))
        
    patches2, texts2 = ax.pie(x = df_counts_opposite['all_met'], labels =  
                None, labeldistance = 0.2, 
                colors = inner_colors, startangle = 180, shadow = False,
                radius=1-size, wedgeprops=dict(width=size, edgecolor='w'))

    labels_sig_group = list(map(lambda x : x[0] + " : " +str(x[1]), zip(sig_group.index, sig_group['counts'])))
    legend1 = plt.legend(patches, labels_sig_group, loc='lower left', prop={'size': 5},  framealpha = 0.3)
    legend2 = plt.legend(patches2, df_counts_opposite['all_met'], loc='upper left', prop={'size': 6}, framealpha = 0.5)
    ax.add_artist(legend1)
    ax.add_artist(legend2)


def myplot(file, xlabel, ylabel, types):
    fig1 = plt.figure(figsize=(8,8), facecolor='white',edgecolor='black')
    grid = plt.GridSpec(7, 7, wspace=0, hspace=0)
 
    ax_x_1 = fig1.add_subplot(grid[1, 0:5])
    ax_x_2 = fig1.add_subplot(grid[0, 0:5])
    ax_y_1 = fig1.add_subplot(grid[2: 7,5])
    ax_y_2 = fig1.add_subplot(grid[2: 7,6])
    
    if types == 1:
        sns.histplot(ax = ax_x_1, x = x, kde=True, stat = 'count', edgecolor = '#E0E0E0')
        plot_scatter(ax_x_2, x, y, data, colors, labels, xlabel, ylabel, 0.05, 0.05, xlim = (0, 1), ylim = (0, 0.05), pointSize=[2, 2], pshape = ['o', 'o'], types = types)
        sns.histplot(ax = ax_y_1, y = y, kde=True, stat = 'count', edgecolor = '#E0E0E0')
        plot_scatter(ax_y_2, x, y, data, colors, labels, xlabel, ylabel, 0.05, 0.05, xlim = (0, 0.05), ylim = (0, 1), pointSize=[2, 2], pshape = ['o', 'o'],  types = types)
    
    elif types == 2:
        sns.histplot(ax = ax_x_1, x = x_same, kde=True, stat = 'count', edgecolor = '#E0E0E0', color = colors[0])	
        sns.histplot(ax = ax_x_2, x = x_oppo, kde=True, stat = 'count', edgecolor = '#E0E0E0', color = colors[1])
        sns.histplot(ax = ax_y_1, y = y_same, kde=True, stat = 'count', edgecolor = '#E0E0E0', color = colors[0])	
        sns.histplot(ax = ax_y_2, y = y_oppo, kde=True, stat = 'count', edgecolor = '#E0E0E0', color = colors[1])

    elif types == 3:
        left_show = np.percentile(x, 10)
        right_show = np.percentile(x, 80)
        up_show = np.percentile(y[~y.isnull()], 0.01)
        bot_show = np.percentile(y[~y.isnull()], 99.99)
        
        bins_x = range(round(left_show), round(right_show), round((abs(left_show)+abs(right_show))/20))
        bins_y = range(round(up_show), round(bot_show), round((abs(up_show)+abs(bot_show))/20))
        sns.histplot(ax = ax_x_1, x = x_same, kde=False, stat = 'density', edgecolor = '#E0E0E0', color = colors[0], bins = bins_x)	
        sns.histplot(ax = ax_x_2, x = x_oppo, kde=False, stat = 'density', edgecolor = '#E0E0E0', color = colors[1], bins = bins_x)
        sns.histplot(ax = ax_y_1, y = y_same, kde=False, stat = 'density', edgecolor = '#E0E0E0', color = colors[0], bins = bins_y)	
        sns.histplot(ax = ax_y_2, y = y_oppo, kde=False, stat = 'density', edgecolor = '#E0E0E0', color = colors[1], bins = bins_y)

    else:
    	print("invalid plot types")
    	return
    
    ax_x_1.get_xaxis().set_visible(False)
    ax_x_1.get_yaxis().set_visible(False)
    ax_y_1.get_xaxis().set_visible(False)
    ax_y_1.get_yaxis().set_visible(False)
    ax_x_2.get_xaxis().set_visible(False)
    ax_y_2.get_yaxis().set_visible(False)
    if types == 2 or types == 3:
        ax_x_2.get_yaxis().set_visible(False)
        ax_y_2.get_xaxis().set_visible(False)
        pass
    
    ax1 = fig1.add_subplot(grid[2:7,0:5])
    if types !=3:
        ax1 = plot_scatter(ax1, x, y, data, colors, labels, xlabel, ylabel, 0.05, 0.05, xlim = (0,1), ylim=(0,1),rec=True, legend=True, pointSize = [1.2, 1.2], pshape = pshape,  types = types)
        ax1.annotate('0.05 shreshold', xy = (0.05, 0.05), xytext= (0.15, 0.2), 
                     xycoords='data', fontsize=9, ha='center', va='center',
                     bbox=dict(boxstyle='square', fc='firebrick'),
                     arrowprops=dict(arrowstyle='->, head_length=0.4,head_width=0.2', lw=2.0, color='black'), 
                     color='white')
    
    else:
        ax1 = plot_scatter(ax1, x, y, data, colors, labels, xlabel, ylabel, 0, 0,
                  xlim = (left_show, right_show),
                  ylim = (up_show, bot_show),
                  rec=False, legend=True, pointSize = [1.2, 1.2, 2, 2, 2], pshape = pshape,  types = types)

    ax_stat_all = fig1.add_subplot(grid[0:2, 5:7])
    sigGroup_color = ['#FF0000', '#00FF00', '#9370DB', '#A9A9A9']
    plot_pie(ax_stat_all, sig_group, df_counts_opposite, sigGroup_color, colors)

    fig1.set_rasterized(True)
    fig1.savefig(file, dpi = 300) 


myplot(file, xlabel, ylabel, types = types)

