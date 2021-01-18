#!/usr/bin/env python
# coding: utf-8

import pandas as pd
import numpy as np
import argparse
import matplotlib.pyplot as plt
import summarize_baseline2methods_classRatio as classratio



def single_read_genes_Stat(dat_sr):
    n = dat_sr.shape[0]
    dat_sr_flag_in_FLAIR = dat_sr[dat_sr["flag_in_FLAIR"] == 1]
    dat_sr_flag_in_isoseq = dat_sr[dat_sr["flag_in_isoseq3Cluster"] == 1]

    removed_by_isoseq3Cluster = n - dat_sr['flag_in_isoseq3Cluster'].sum()
    removed_by_flair = n - dat_sr['flag_in_FLAIR'].sum()
    removed_by_sqantiQC = n - dat_sr['flag_in_SQANTI3filter'].sum()

    Internal_modified_by_isoseq3Cluster = dat_sr_flag_in_isoseq['flag_modified_internal_in_isoseq3Cluster'].sum()
    Internal_modified_by_flair = dat_sr_flag_in_FLAIR['flag_modified_internal_in_FLAIR'].sum()
    Internal_modified_by_SQANTIfilter = 0

    single_read_cluster_isoseq3Cluster = dat_sr_flag_in_isoseq['flag_single_read_cluster_isoseq3Cluster'].sum()
    single_read_cluster_FLAIR = dat_sr_flag_in_FLAIR['flag_single_read_cluster_FLAIR'].sum()
    single_read_cluster_SQANTI = None
    
    res = pd.DataFrame({'flair': [n, removed_by_flair, n - removed_by_flair ,Internal_modified_by_flair, single_read_cluster_FLAIR], 
                    'isoseq3Cluster': [n, removed_by_isoseq3Cluster, n - removed_by_isoseq3Cluster, Internal_modified_by_isoseq3Cluster, single_read_cluster_isoseq3Cluster], 
                    'SQANTIfilter': [n, removed_by_sqantiQC , n - removed_by_sqantiQC, Internal_modified_by_SQANTIfilter, single_read_cluster_SQANTI]})
    res.index = ['No. Total', 'No. Removed', 'No. Retained','No. Internal_modified (within retained)', 'No. Single_read_cluster (within retained)']
    res.index.name = 'single_read_genes'
    return res


def make_df_classCounts(data, name):
    if name == 'total':
        data_in_name = data
    else:
        col = 'flag_in_' + name
        data_in_name = data[data[col] == 0]
    data_clssCounts = data_in_name.groupby('structural_category')["readID"].count()
    data_clssCounts.name = name
    data_clssCounts = pd.DataFrame(data_clssCounts)
    return data_clssCounts


def single_read_gene_cls_Stat(dat_sr):
    list_sr_flag_notin_isoseq_cls = make_df_classCounts(dat_sr, 'isoseq3Cluster')
    list_sr_flag_notin_FLAIR_cls = make_df_classCounts(dat_sr, 'FLAIR')
    list_sr_flag_notin_SQANTI_cls = make_df_classCounts(dat_sr, 'SQANTI3filter')
    dat_sr_total_cls = make_df_classCounts(dat_sr, 'total')
    
    df_sr_removeStat = dat_sr_total_cls.merge(
    list_sr_flag_notin_isoseq_cls, left_index = True, right_index = True, how = 'outer').merge(
    list_sr_flag_notin_FLAIR_cls, left_index = True, right_index = True, how = 'outer').merge(
    list_sr_flag_notin_SQANTI_cls, left_index = True, right_index = True, how = 'outer')
    
    df_sr_removeStat['isoseq3Cluster_ratio'] = df_sr_removeStat['isoseq3Cluster']/df_sr_removeStat['total']
    df_sr_removeStat['FLAIR_ratio'] = df_sr_removeStat['FLAIR']/df_sr_removeStat['total']
    df_sr_removeStat['SQANTI3filter_ratio'] = df_sr_removeStat['SQANTI3filter']/df_sr_removeStat['total']
    df_sr_removeStat = df_sr_removeStat[['total', 'isoseq3Cluster', 'isoseq3Cluster_ratio', 'FLAIR', 'FLAIR_ratio', 'SQANTI3filter', 'SQANTI3filter_ratio']]
    if not 'Accumulative' in df_sr_removeStat.index:
        df_sr_removeStat.loc['Accumulative'] = df_sr_removeStat.sum(axis=0)
        for i in df_sr_removeStat.columns:
            if 'ratio' in i:
                df_sr_removeStat.loc['Accumulative'][i] = None
    return df_sr_removeStat


def make_geneCluster_numTable(dat):
    iso_cluster = pd.DataFrame(dat[(dat['flag_multigene_isoseq3Cluster'] == 0) & (dat['flag_in_isoseq3Cluster'] == 1)].groupby('gene_id')['iso_cluster_id'].nunique())
    iso_cluster.columns = ['No. cluster']
    flair_cluster = pd.DataFrame(dat[(dat['flag_multigene_FLAIR']==0) & (dat['flag_in_FLAIR'] == 1)].groupby('gene_id')['flair_transcript_id'].nunique())
    flair_cluster.columns = ['No. cluster']
    clusterStat = iso_cluster.merge(flair_cluster, left_index=True, right_index=True, how = 'outer', suffixes = ('_iso', '_flair'))
    clusterStat['flag_isoGreater'] = np.where(clusterStat['No. cluster_iso'] > clusterStat['No. cluster_flair'], 1, 0) 
    clusterStat['flag_flairGreater'] = np.where(clusterStat['No. cluster_iso'] < clusterStat['No. cluster_flair'], 1, 0)
    clusterStat['flag_equalNum'] = np.where(clusterStat['No. cluster_iso'] == clusterStat['No. cluster_flair'], 1, 0)
    clusterStat['flag_ClusterInFLAIR_notInIsoseq'] = np.where(clusterStat['No. cluster_iso'].isna(), 1, 0)
    clusterStat['flag_ClusterInIsoseq_notInFLAIR'] = np.where(clusterStat['No. cluster_flair'].isna(), 1, 0)
    return clusterStat


def mthan2read_classCounts(dat, method, criteria, check_fsm = None, equal_list = pd.DataFrame(), flag_in = True):
    data = dat[dat['flag_single_read_gene_baseline'] == 0]
    flag_multigene = 'flag_multigene_' + method
    flag_in_method = 'flag_in_' + method
    flag_if_modify = 'flag_modified_internal_in_' + method
    if criteria == 'EqualOneCluster' or criteria == 'EqualMoreThanOne':
        assert not equal_list.empty  
    
    if method == 'isoseq3Cluster':
        idcol = 'iso_cluster_id'
    elif method == 'FLAIR':
        idcol = 'flair_transcript_id'  
    
    if flag_in:
        dat_mt2read = data[(data[flag_multigene] == 0) & (data[flag_in_method] == 1)]
    else:
        dat_mt2read = data[data[flag_multigene] == 0]

    if criteria == 'OneCluster':
        dat_mt2read_num = dat_mt2read.groupby('gene_id')[idcol].nunique()
        dat_mt2read_geneList = dat_mt2read_num[dat_mt2read_num == 1].index
    elif criteria == 'MoreThanOne':
        dat_mt2read_num = dat_mt2read.groupby('gene_id')[idcol].nunique()
        dat_mt2read_geneList = dat_mt2read_num[dat_mt2read_num > 1].index
    elif criteria == 'EqualOneCluster' or criteria == 'EqualMoreThanOne':
        dat_mt2read_geneList = equal_list
    
    dat_mt2read_clspg = dat_mt2read[dat_mt2read['gene_id'].isin(dat_mt2read_geneList)]
        
    if check_fsm == 'AtLeastOneEither':
        dat_mt2read_geneList = dat_mt2read_clspg[dat_mt2read_clspg['structural_category'] == 'full-splice_match'].drop_duplicates('gene_id')['gene_id']
        dat_mt2read_clspg = dat_mt2read_clspg[dat_mt2read_clspg['gene_id'].isin(dat_mt2read_geneList)]
    elif check_fsm == 'AtLeastOneBoth':
        DictCLSwithFSM = dict(dat_mt2read_clspg.groupby(idcol)['structural_category'].unique().apply(
            lambda x: 1 if 'full-splice_match' in x else 0).items())
        temp = dat_mt2read_clspg.groupby('gene_id')[idcol].unique().apply(
            check_cluster_status, args = (DictCLSwithFSM,))
        dat_mt2read_geneList = temp[temp==1].index
        dat_mt2read_clspg = dat_mt2read_clspg[dat_mt2read_clspg['gene_id'].isin(dat_mt2read_geneList)]
    
    #modified, notmodified = make_dfs_isModified_byIfmodify(dat_mt2read_clspg, method)
    return dat_mt2read_clspg


def check_cluster_status(clusterIDs, clusterDict):
    status = 1
    if len(clusterIDs) == 0:
        return 0
    for i in clusterIDs:
        if (i in clusterIDs) & (i != None):
            if clusterDict[i] == 0:
                status = 0
                break
    return status


def make_dfs_isModified_byIfmodify(data, method):
    df = data
    modifyFlag = 'flag_modified_internal_in_' + method
    modified = df[df[modifyFlag] == 1].groupby('structural_category')['readID'].count()
    notmodified = df[df[modifyFlag] == 0].groupby('structural_category')['readID'].count()
    modified.index.name = 'type'
    modified.name = 'modified_by_' + method.replace('isoseq3Cluster', 'isoseq3')
    notmodified.index.name = 'type'
    notmodified.name = 'not_modified_by_' + method.replace('isoseq3Cluster', 'isoseq3')
    return pd.DataFrame(modified), pd.DataFrame(notmodified)


def make_dfs_isModified_byClass(data):
    df = data
    df_modified_both = pd.concat([df[(df['flag_modified_internal_in_FLAIR'] == 1) & 
                                      (df['flag_in_FLAIR'] == 1) & 
                                      (df['flag_modified_internal_in_isoseq3Cluster'] == 1) & 
                                      (df['flag_in_isoseq3Cluster'] == 1)],    df[(df['flag_in_FLAIR'] == 0) & 
                                                                               (df['flag_in_isoseq3Cluster'] == 0)]])
    df_notmodified_both = df[(df['flag_modified_internal_in_FLAIR'] == 0) & 
                             (df['flag_in_FLAIR'] == 1) & 
                             (df['flag_modified_internal_in_isoseq3Cluster'] == 0) & 
                             (df['flag_in_isoseq3Cluster'] == 1)]
    
    df_modified_isoOnly = pd.concat([df[(df['flag_modified_internal_in_isoseq3Cluster'] == 1) & 
                                      (df['flag_in_isoseq3Cluster'] == 1) & 
                                      (df['flag_modified_internal_in_FLAIR'] == 0) & 
                                      (df['flag_in_FLAIR'] == 1)],    df[(df['flag_in_FLAIR'] == 1) & 
                                                                      (df['flag_in_isoseq3Cluster'] == 0)]])
    
    df_modified_flairOnly = pd.concat([df[(df['flag_modified_internal_in_isoseq3Cluster'] == 0) & 
                                      (df['flag_in_isoseq3Cluster'] == 1) & 
                                      (df['flag_modified_internal_in_FLAIR'] == 1) & 
                                      (df['flag_in_FLAIR'] == 1)],    df[(df['flag_in_FLAIR'] == 0) & 
                                                                      (df['flag_in_isoseq3Cluster'] == 1)]])
    
    m1 = df_modified_both.groupby('structural_category')['readID'].count()
    m2 = df_notmodified_both.groupby('structural_category')['readID'].count()
    m3 = df_modified_isoOnly.groupby('structural_category')['readID'].count()
    m4 = df_modified_flairOnly.groupby('structural_category')['readID'].count()
    m1.name = 'modified_both'
    m2.name = 'not_modified_both'
    m3.name = 'modified_isoOnly'
    m4.name = 'modified_flairOnly'
    return pd.DataFrame(m1), pd.DataFrame(m2), pd.DataFrame(m3), pd.DataFrame(m4)


def merge_Counts(data, criteria = 'OneCluster', check_FSM = None, groupby = 'ifmodify', clusterStat = pd.DataFrame()):
    if groupby == 'class':
        flag_in = False
    else:
        flag_in = True
        
    if criteria == 'OneCluster' or criteria == 'MoreThanOne':
        subdata_matchCriteria_flair = mthan2read_classCounts(data, 'FLAIR', criteria = criteria, check_fsm = check_FSM,  flag_in = flag_in)
        subdata_matchCriteria_iso = mthan2read_classCounts(data, 'isoseq3Cluster', criteria = criteria, check_fsm = check_FSM, flag_in = flag_in)
    elif criteria == 'EqualOneCluster' or criteria == 'EqualMoreThanOne':
        assert not clusterStat.empty
        if criteria == 'EqualOneCluster':
            condition = clusterStat[clusterStat['flag_equalNum'] == 1]['No. cluster_iso'] == 1
        elif criteria == 'EqualMoreThanOne':
            condition = clusterStat[clusterStat['flag_equalNum'] == 1]['No. cluster_iso'] > 1
        equal_list_passin = clusterStat[clusterStat['flag_equalNum'] == 1][condition].index
        
        subdata_matchCriteria_flair = mthan2read_classCounts(data, 'FLAIR', criteria = criteria, check_fsm = check_FSM, equal_list = equal_list_passin, flag_in = flag_in)
        subdata_matchCriteria_iso = mthan2read_classCounts(data, 'isoseq3Cluster', criteria = criteria, check_fsm = check_FSM, equal_list = equal_list_passin, flag_in = flag_in)      
        
    if not groupby == 'class':
            m1, m2 = make_dfs_isModified_byIfmodify(subdata_matchCriteria_flair, 'FLAIR')
            m3, m4 = make_dfs_isModified_byIfmodify(subdata_matchCriteria_iso, 'isoseq3Cluster')
    else:
        subdata_matchCriteria_both = pd.concat([subdata_matchCriteria_flair, subdata_matchCriteria_iso]).drop_duplicates()
        m1, m2, m3, m4 = make_dfs_isModified_byClass(subdata_matchCriteria_both)
        
    m = m1
    for i in [m2, m3, m4]:
        m = m.merge(i, left_index = True, right_index = True, how = 'outer')
    m = m.fillna(0)
    if groupby == 'class':
        m = m.T
    return m, classratio.generate_dfRatio(m)


def df_withN_supportingReads(data, shreshold = 2, largerThan = False):
    df = data
    gene_counts = df.groupby('gene_id')['readID'].count()
    if not largerThan:
        geneMatchList = gene_counts[gene_counts == shreshold].index
    else:
        geneMatchList = gene_counts[gene_counts > shreshold].index
    res = data[data['gene_id'].isin(geneMatchList)]
    return res


def make_df_clusterCompare_BysupReads(rangefrom, rangeto, data):
    flag_isoGreater = []
    flag_flairGreater = []
    flag_equalNum = []
    flag_ClusterInFLAIR_notInIsoseq = []
    flag_ClusterInIsoseq_notInFLAIR = []
    category = ['flag_isoGreater', 
                'flag_flairGreater', 
                'flag_equalNum','flag_ClusterInFLAIR_notInIsoseq', 
                'flag_ClusterInIsoseq_notInFLAIR']
    categoryList = [flag_isoGreater, 
                    flag_flairGreater, 
                    flag_equalNum, 
                    flag_ClusterInFLAIR_notInIsoseq, 
                    flag_ClusterInIsoseq_notInFLAIR]
    
    for i in range(rangefrom, rangeto+1):
        df = df_withN_supportingReads(data, i)
        clusterStat_table = make_geneCluster_numTable(df)
        clusterStat_out = clusterStat_table[[i for i in clusterStat_table.columns if 'flag' in i]].sum().to_frame()
        clusterStat_out.columns = ['counts']
        #outfile_clusterCount = outdir + '_geneWith_' + str(i) + '_supportReads_ClusterPerGene_counts.txt'
        s = 0
        for j in category:
            categoryList[s].append(int(clusterStat_out.loc[j]))
            s+=1

        if i == rangeto:
            df = df_withN_supportingReads(data, i, largerThan = True)
            clusterStat_table = make_geneCluster_numTable(df)
            clusterStat_out = clusterStat_table[[i for i in clusterStat_table.columns if 'flag' in i]].sum().to_frame()
            clusterStat_out.columns = ['counts']
            s = 0
            for j in category:
                categoryList[s].append(int(clusterStat_out.loc[j]))
                s+=1

    df_sum = pd.DataFrame({'flag_isoGreater' : flag_isoGreater, 
              'flag_flairGreater' : flag_flairGreater, 
              'flag_equalNum' : flag_equalNum, 
              'flag_ClusterInFLAIR_notInIsoseq' : flag_ClusterInFLAIR_notInIsoseq, 
              'flag_ClusterInIsoseq_notInFLAIR': flag_ClusterInIsoseq_notInFLAIR})
    df_sum.index = [i for i in range(rangefrom, rangeto + 1)] + ['Greater than ' + str(rangeto)]
    df_sum.index.name = 'No. supporting read'
    return df_sum




