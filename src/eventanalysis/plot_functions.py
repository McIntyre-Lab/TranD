#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 26 11:33:04 2021

@author: adalena.nanni
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
import seaborn as sns
from upsetplot import UpSet

## Plot Functions
   
def plot_transcript_in_gene_pie(md_data,f1_odds,f2_odds,name1,name2):
    # Get counts
    geneCount = md_data[['gene_id','transcript_in_gene']].drop_duplicates()['transcript_in_gene'].value_counts(sort=False)
    geneCount[name1+'_only'] = f1_odds['gene_id'].nunique()
    geneCount[name2+'_only'] = f2_odds['gene_id'].nunique()
    # Make dictionary of names for plot and order
    geneCountNameDict = {'equal':'Equal',name1+'_greater':name1+" Greater",name2+'_greater':name2+" Greater",name1+'_only':name1+" Only",name2+'_only':name2+" Only"}
    geneCountOrderDict = {'equal':0,name1+'_greater':1,name2+'_greater':2,name1+'_only':3,name2+'_only':4}
    # Reindex to order properly
    geneCount = geneCount.reindex(list(geneCountOrderDict))
    # Plot pie chart
    geneCount.plot(kind='pie',y='transcript_in_gene',labels=None,figsize=(12,6),autopct=(lambda pct: get_pie_label(pct,geneCount)),colormap=ListedColormap(sns.color_palette('colorblind',15).as_hex()))
    plt.legend(bbox_to_anchor=(0.85,1), loc="upper left",labels=geneCount.index.map(geneCountNameDict))
    plt.ylabel("")
    plt.tight_layout()

def plot_transcript_in_gene_split_pie(md_data,f1_odds,f2_odds,name1,name2):
    # Get counts
    geneF1only = f1_odds.groupby('gene_id')['transcript_id'].nunique().reset_index().rename(columns={'transcript_id':'num_transcript_in_gene_'+name1})
    geneF1only['num_transcript_in_gene_'+name2] = 0
    geneF1only['transcript_in_gene'] = name1+'_only'
    geneF2only = f2_odds.groupby('gene_id')['transcript_id'].nunique().reset_index().rename(columns={'transcript_id':'num_transcript_in_gene_'+name2})
    geneF2only['num_transcript_in_gene_'+name1] = 0
    geneF2only['transcript_in_gene'] = name2+'_only'

    # Make dictionary of names for plot and order
    geneCountNameDict = {'equal':'Equal',name1+'_greater':name1+" Greater",name2+'_greater':name2+" Greater",name1+'_only':name1+" Only",name2+'_only':name2+" Only"}
    geneCountOrderDict = {'equal':0,name1+'_greater':1,name2+'_greater':2,name1+'_only':3,name2+'_only':4}

    # Plot pie chart
    # 1) all genes
    # 2) genes with 1 xcrpt for at least 1 dataset
    # 3) genes with >1 xcrpt for at least one dataset
    # 4) genes with > 1 xcrpt for both
    fig = plt.figure()
    for num in range(1,5):
        ax = plt.subplot2grid((2,2),(abs(num%2-1),int(num/3)),fig=fig)
        if num == 1:
            geneCount = pd.concat([md_data[['gene_id','transcript_in_gene']].drop_duplicates(),
                                   geneF1only,geneF2only],ignore_index=True,sort=False)['transcript_in_gene'].value_counts(sort=False)
            title = "All Genes (n = {})".format(geneCount.sum())
            legendLabels = geneCount.reindex(list(geneCountOrderDict)).index.map(geneCountNameDict)
        elif num == 2:
            geneCount = pd.concat([md_data[(md_data['num_transcript_in_gene_'+name1]==1)|(md_data['num_transcript_in_gene_'+name2]==1)][['gene_id','transcript_in_gene']].drop_duplicates(),
                                   geneF1only[geneF1only['num_transcript_in_gene_'+name1]==1],
                                   geneF2only[geneF2only['num_transcript_in_gene_'+name2]==1]],ignore_index=True,sort=False)['transcript_in_gene'].value_counts(sort=False)
            title = "Genes with 1 Transcript\n(in at least one dataset, n = {})".format(geneCount.sum())
        elif num == 3:
            geneCount = pd.concat([geneF1only[geneF1only['num_transcript_in_gene_'+name1]>1],
                                   geneF2only[geneF2only['num_transcript_in_gene_'+name2]>1]],ignore_index=True,sort=False)['transcript_in_gene'].value_counts(sort=False)
            title = "Genes with >1 Transcript\n(exclusive to one dataset, n = {})".format(geneCount.sum())
        else:
            geneCount = md_data[(md_data['num_transcript_in_gene_'+name1]>1)&(md_data['num_transcript_in_gene_'+name2]>1)][['gene_id','transcript_in_gene']].drop_duplicates()['transcript_in_gene'].value_counts(sort=False)
            title = "Genes with >1 Transcript\n(in both datasets, n = {})".format(geneCount.sum())            
        geneCount.reindex(list(geneCountOrderDict)).plot(kind='pie',y='transcript_in_gene',labels=None,figsize=(12,6),autopct=(lambda pct: get_pie_label(pct,geneCount)),colormap=ListedColormap(sns.color_palette('colorblind',15).as_hex()),ax=ax)
        ax.set_ylabel("")
        ax.set_title(title)
    plt.legend(bbox_to_anchor=(1.2,1.2), loc="upper left",labels=legendLabels)
    plt.tight_layout()

def plot_gene_stack(md_data,name1,name2,useProp=False):
    genePairDF = md_data.groupby(['transcript_in_gene','recip_min_pair_in_gene'])['gene_id'].nunique().reset_index().pivot_table(index=['transcript_in_gene'],columns='recip_min_pair_in_gene',values='gene_id')
    genePairDF = genePairDF.rename(columns={'match':'Match',
                                            'partial_match':'Partial Match',
                                            'no_match':'No Match',
                                            'subset':'Subset',
                                            'partial_subset':'Partial Subset',
                                            'no_subset':'No Subset'},
                                    index={name1+'_greater':name1+" Greater",
                                           name2+'_greater':name2+" Greater",
                                           "equal":"Equal"})
    if useProp:
        genePairDF = genePairDF.div(genePairDF.sum(axis=1), axis=0)
        title = "Proportion of Genes"
    else:
        title = "Number of Genes"
    genePairDF.plot(kind='bar',figsize=(10,6),stacked=True,rot=45,colormap=ListedColormap(sns.color_palette('Paired').as_hex()))
    plt.ylabel(title)
    plt.xlabel("")
    plt.legend(bbox_to_anchor=(1.04,1), loc="upper left")
    plt.tight_layout()

def get_pie_label(pct, allvals):
    absolute = int(round(pct/100.*np.sum(allvals)))
    return "{:.1f}%\n({:d})".format(pct, absolute)

def plot_transcript_in_gene_upset(md_data,f1_odds,f2_odds,name1,name2):
    geneF1only = f1_odds.groupby('gene_id')['transcript_id'].nunique().reset_index().rename(columns={'transcript_id':'num_transcript_in_gene_'+name1})
    geneF1only['num_transcript_in_gene_'+name2] = 0
    geneF1only['transcript_in_gene'] = name1+'_only'
    geneF2only = f2_odds.groupby('gene_id')['transcript_id'].nunique().reset_index().rename(columns={'transcript_id':'num_transcript_in_gene_'+name2})
    geneF2only['num_transcript_in_gene_'+name1] = 0
    geneF2only['transcript_in_gene'] = name2+'_only'
    geneAll = pd.concat([md_data[['gene_id','num_transcript_in_gene_'+name1,'num_transcript_in_gene_'+name2,'transcript_in_gene']].drop_duplicates(),
                         geneF1only,geneF2only],ignore_index=True,sort=False)
    colList = []
    # Flag when genes have 1-4 transcripts
    for num in range(1,5):
        geneAll[str(num)+" Transcript(s) "+name1] = np.where(geneAll['num_transcript_in_gene_'+name1]==num,True,False)
        colList.append(str(num)+" Transcript(s) "+name1)
        geneAll[str(num)+" Transcript(s) "+name2] = np.where(geneAll['num_transcript_in_gene_'+name2]==num,True,False)
        colList.append(str(num)+" Transcript(s) "+name2)
    # Flag when genes have 5+ transcripts
    geneAll["5+ Transcript(s) "+name1] = np.where(geneAll['num_transcript_in_gene_'+name1]>=5,True,False)
    colList.append("5+ Transcript(s) "+name1)
    geneAll["5+ Transcript(s) "+name2] = np.where(geneAll['num_transcript_in_gene_'+name2]>=5,True,False)
    colList.append("5+ Transcript(s) "+name2)
    plot_upset(geneAll.set_index(colList),"Number of Genes with Each Number of Transcripts")
        
def plot_recip_min_pair_AS_upset(md_data):
    recipMinPairAS = md_data[md_data['flag_recip_min_match']==1][['transcript_1','transcript_2',
                            'flag_alt_exon_recip_min_match','flag_alt_donor_acceptor_recip_min_match',
                            'flag_IR_recip_min_match','flag_5_variation_recip_min_match',
                            'flag_3_variation_recip_min_match','flag_nonoverlapping_recip_min_match']].copy()
    recipMinPairAS = recipMinPairAS.rename(columns={'flag_alt_exon_recip_min_match':'Alt. Exon',
                                                    'flag_alt_donor_acceptor_recip_min_match':'Alt. Donor/Acceptor',
                                                    'flag_IR_recip_min_match':'Intron Retention',
                                                    'flag_5_variation_recip_min_match':'5\' Variation',
                                                    'flag_3_variation_recip_min_match':'3\' Variation',
                                                    'flag_nonoverlapping_recip_min_match':'Nonoverlapping'})
    plot_upset(recipMinPairAS.set_index([c for c in recipMinPairAS.columns if "transcript" not in c]),
               "Number of Reciprocal Minimum Pairs with AS Categories")

def plot_recip_min_pair_AS_upset_nt_box(md_data):
    recipMinPairAS = md_data[md_data['flag_recip_min_match']==1][['transcript_1','transcript_2',
                            'flag_alt_exon_recip_min_match','flag_alt_donor_acceptor_recip_min_match',
                            'flag_IR_recip_min_match','flag_5_variation_recip_min_match',
                            'flag_3_variation_recip_min_match','flag_nonoverlapping_recip_min_match',
                            'num_nt_diff','prop_nt_diff']].copy()
    recipMinPairAS = recipMinPairAS.rename(columns={'flag_alt_exon_recip_min_match':'Alt. Exon',
                                                    'flag_alt_donor_acceptor_recip_min_match':'Alt. Donor/Acceptor',
                                                    'flag_IR_recip_min_match':'Intron Retention',
                                                    'flag_5_variation_recip_min_match':'5\' Variation',
                                                    'flag_3_variation_recip_min_match':'3\' Variation',
                                                    'flag_nonoverlapping_recip_min_match':'Nonoverlapping',
                                                    'num_nt_diff':'Number NT Different',
                                                    'prop_nt_diff':'Proportion NT Different'})
    plot_upset(recipMinPairAS.set_index([c for c in recipMinPairAS.columns if ("transcript" not in c) and ("NT" not in c)]),
               "Number of Reciprocal Minimum Pairs with AS Categories",['Number NT Different','Proportion NT Different'])

def plot_gene_AS_upset(md_data):
    geneRecipMatchAS = md_data.groupby('gene_id')[['flag_alt_exon_recip_min_match','flag_alt_donor_acceptor_recip_min_match',
                                                   'flag_IR_recip_min_match','flag_5_variation_recip_min_match',
                                                   'flag_3_variation_recip_min_match',
                                                   'flag_nonoverlapping_recip_min_match']].max().astype(bool).reset_index()
    geneRecipMatchAS = geneRecipMatchAS.rename(columns={'flag_alt_exon_recip_min_match':'Alt. Exon',
                                                        'flag_alt_donor_acceptor_recip_min_match':'Alt. Donor/Acceptor',
                                                        'flag_IR_recip_min_match':'Intron Retention',
                                                        'flag_5_variation_recip_min_match':'5\' Variation',
                                                        'flag_3_variation_recip_min_match':'3\' Variation',
                                                        'flag_nonoverlapping_recip_min_match':'Nonoverlapping'})
    plot_upset(geneRecipMatchAS.set_index([c for c in geneRecipMatchAS.columns if c != "gene_id"]),
               "Number of Genes with AS Categories in Reciprocal Minimum Pairs")

def plot_gene_AS_upset_nt_box(md_data,name1,name2):
    recipMinPairAS = md_data[md_data['flag_recip_min_match']==1][['gene_id',
                            'flag_alt_exon_recip_min_match','flag_alt_donor_acceptor_recip_min_match',
                            'flag_IR_recip_min_match','flag_5_variation_recip_min_match',
                            'flag_3_variation_recip_min_match','flag_nonoverlapping_recip_min_match',
                            'num_recip_min_match_in_gene','num_transcript_in_gene_'+name1,
                            'num_transcript_in_gene_'+name2,'num_nt_diff','prop_nt_diff']].copy()
    geneRecipMatchAS = recipMinPairAS.groupby('gene_id').agg({'flag_alt_exon_recip_min_match':'max',
                                                              'flag_alt_donor_acceptor_recip_min_match':'max',
                                                              'flag_IR_recip_min_match':'max',
                                                              'flag_5_variation_recip_min_match':'max',
                                                              'flag_3_variation_recip_min_match':'max',
                                                              'flag_nonoverlapping_recip_min_match':'max',
                                                              'num_recip_min_match_in_gene':'max',
                                                              'num_transcript_in_gene_'+name1:'max',
                                                              'num_transcript_in_gene_'+name2:'max',
                                                              'num_nt_diff':'mean',
                                                              'prop_nt_diff':'mean'}).reset_index()
    geneFlagCols = ['flag_alt_exon_recip_min_match','flag_alt_donor_acceptor_recip_min_match',
                    'flag_IR_recip_min_match','flag_5_variation_recip_min_match',
                    'flag_3_variation_recip_min_match','flag_nonoverlapping_recip_min_match']
    geneRecipMatchAS[geneFlagCols] = geneRecipMatchAS[geneFlagCols].astype(bool)
    geneRecipMatchAS = geneRecipMatchAS.rename(columns={'flag_alt_exon_recip_min_match':'Alt. Exon',
                                                        'flag_alt_donor_acceptor_recip_min_match':'Alt. Donor/Acceptor',
                                                        'flag_IR_recip_min_match':'Intron Retention',
                                                        'flag_5_variation_recip_min_match':'5\' Variation',
                                                        'flag_3_variation_recip_min_match':'3\' Variation',
                                                        'flag_nonoverlapping_recip_min_match':'Nonoverlapping',
                                                        'num_recip_min_match_in_gene':'# Recip. Min.\nMatch Transcripts',
                                                        'num_transcript_in_gene_'+name1:'# Transcripts\nin '+name1,
                                                        'num_transcript_in_gene_'+name2:'# Transcripts\nin '+name2,
                                                        'num_nt_diff':'Avg #\nNT Different',
                                                        'prop_nt_diff':'Avg Proportion\nNT Different'})
    plot_upset(geneRecipMatchAS.set_index([c for c in geneRecipMatchAS.columns if ("gene" not in c) and ("NT" not in c) and ("Transcripts" not in c) ]),
               "",['# Recip. Min.\nMatch Transcripts','# Transcripts\nin '+name1,'# Transcripts\nin '+name2,'Avg #\nNT Different','Avg Proportion\nNT Different'])

def plot_upset(df,title,boxCols=None):
    upset = UpSet(df,subset_size='count',show_counts=True,sort_by='cardinality',sort_categories_by='cardinality')
    if boxCols is not None:
        if type(boxCols) != list:
            boxCols = list(boxCols)
        for col in boxCols:
            upset.add_catplot(value=col, kind='box', elements = 4, showfliers=False, color=sns.color_palette('colorblind',15).as_hex()[boxCols.index(col)])
    upset.plot()
    plt.subplots_adjust(right=1.00001)
    plt.suptitle(title)
