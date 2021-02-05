#!/usr/bin/env python

# Functions for FSM Consolidation

import pandas as pd

def split_transcript_id(df,col_name=None,sort_list=None):
    # Split transcript_id by pipes and keep all other values the same
    # Some events are found in multiple transcripts (common or constitutive) and
    #     have a piped list of transcript_id values in the transcript_id column
    if col_name == None:
        col_name = 'transcript_id'
    splitList = df[col_name].str.split("|").apply(pd.Series, 1).stack()
    splitList.index = splitList.index.droplevel(-1)
    tempDF = df.copy()
    del(tempDF[col_name])
    splitDF = tempDF.join(splitList.rename(col_name))
    if sort_list != None:
        splitDF = splitDF.sort_values(by=sort_list)
    del(tempDF, splitList)
    return splitDF

def collapse_junctions(df):
    # Make transcript-level junction variables by collapsing by transcript_id
    collapseJuncDF = df.groupby('transcript_id').apply(
            func=lambda x:"|".join(x['junction_id'])).reset_index().rename(
                    columns={0:'junctionID_order'}).sort_values(['junctionID_order'])
    return collapseJuncDF

def collapse_fragments(df, sort_list=None, level=None, featureType=None, featureName=None):
    # Sort values if list of variables to sort by is set
    if sort_list != None:
        df = df.sort_values(by=sort_list)
    # Set feature type input (fragment or fusion) to set id column name
    # Default is "fragment"
    if featureType == None:
        featureType = "fragment"
    if featureName == None:
        featureName = "fragment"
    # Set level of collaping: gene or transcript
    if level == None:
        level = "transcript"
        otherLevel = "gene"
    elif level == "transcript":
        otherLevel = "gene"
    elif level == "gene":
        otherLevel = "transcript"
    else:
        raise Exception('Function "collapse_fragments" requires level definition of \"transcript\" or \"gene\"')
    # Make transcript-level fragment variables by collapsing by transcript_id
    collapseDF = df.groupby(level+'_id').agg({
            otherLevel+"_id":[lambda x:"|".join(x.unique())],
            featureType+'_id':[lambda x:"|".join(x),'count'],
            'chr':['first'],
            featureType+'_start':['min'],
            featureType+'_stop':['max'],
            featureName+'_length':['sum',lambda x:"|".join(x.map(str))]}).reset_index()
    collapseDF.columns = collapseDF.columns.droplevel(1)
    # Rename duplicated column names by index
    collapseDF.columns.values[3] = "num_"+featureName
    collapseDF.columns.values[8] = featureName+"Length_order"
    # Rename columns after aggregation functions
    collapseDF = collapseDF.rename(columns={
            featureType+'_id':featureName+'ID_order',
            featureType+'_start':'start',
            featureType+'_stop':'end',
            featureName+'_length':'transcript_length'})
    if otherLevel == "transcript":
        collapseDF = collapseDF.rename(columns={'transcript_id':'transcriptID_cat'})
    return collapseDF

def modify_gff_attributes(feature, row, diff_start=None, diff_end=None):
    # Modify gff attributes to match the FSM consolidation
    # Change feature source
    feature.source = "FSM_consolidation"
    # Change the transcript ID value to the FSM_consolidation_transcript_id
    if 'ID' in feature.attributes:
        feature.attributes['ID'] = row['FSM_consolidation_transcript_id']
    # If ther feature is an exon or intron, change Parent to be FSM_consolidation_transcript_id
    if feature.featuretype == 'exon' or feature.featuretype == 'intron':
        feature.attributes['Parent'] = row['FSM_consolidation_transcript_id']
        # Retain name for exons/introns that do not change
        # Modify name for exons that have modified start/end
        if 'Name' in feature.attributes:
            if diff_start != None or diff_end != None:
                feature.attributes['Name'][0] = feature.attributes['Name'][0] + "_ex"
    # Set start or end if the feature requires a modified start or end
    if diff_start != None:
        feature.start = diff_start
    if diff_end != None:
        feature.end = diff_end
    return feature

def add_exon(exon, exonDF):
    # Check if exon is already present in another transcript
    # If yes then add new transcript to Parent list
    # If no then add exon to exon list
    exonSeries = pd.Series(str(exon).split(),index=['chr','source',
                           'feature','start','end','score','strand',
                           'phase','attributes'])
    exonSeries = exonSeries.append(pd.Series(exonSeries['attributes'].split(";"),
                                index=['attribute_Name','attribute_Parent',
                                       'attribute_parent_type']))
    exonSeries = exonSeries.drop('attributes')
    if exonSeries['attribute_Name'] in exonDF['attribute_Name'].values:
        oldParent = exonDF.loc[exonDF['attribute_Name']==exonSeries['attribute_Name'],'attribute_Parent'].values[0]
        newParent = oldParent + "," + exonSeries['attribute_Parent'].split("=")[1]
        exonDF.loc[exonDF['attribute_Name']==exonSeries['attribute_Name'],'attribute_Parent'] = newParent
    else:
        exonDF = exonDF.append(exonSeries,ignore_index=True)
    return exonDF
