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
   
def plot_transcript_in_gene_pie(md_data,f1_odds,f2_odds,name1,name2,legendOut):
    # Get counts
    geneCount = md_data[['gene_id','transcript_in_gene']].drop_duplicates()['transcript_in_gene'].value_counts(sort=False)
    geneCount[name1+'_only'] = f1_odds['gene_id'].nunique()
    geneCount[name2+'_only'] = f2_odds['gene_id'].nunique()
    # Make dictionary of names for plot and order
    geneCountNameDict = {'match':'Match',name1+'_greater':name1+" Greater",name2+'_greater':name2+" Greater",name1+'_only':name1+" Only",name2+'_only':name2+" Only"}
    geneCountOrderDict = {'match':0,name1+'_greater':1,name2+'_greater':2,name1+'_only':3,name2+'_only':4}
    # Reindex to order properly
    geneCount = geneCount.reindex(list(geneCountOrderDict))
    # Plot pie chart
    geneCount.plot(kind='pie',y='transcript_in_gene',labels=None,figsize=(12,6),autopct=(lambda pct: get_pie_label(pct,geneCount)),colormap=ListedColormap(sns.color_palette('colorblind',15).as_hex()))
    plt.legend(bbox_to_anchor=(0.85,1), loc="upper left",labels=geneCount.index.map(geneCountNameDict))
    plt.ylabel("")
#    plt.text(0,-1.25,"Percent of genes with equal number of transcripts between {} and {} (blue), more transcripts in {} compared to {} (dark orange), more transcripts in {} compared to {} (gray), and genes exclusive to {} (light orange) or {} (pink). {} genes total.".format(
#            name1,name2,name1,name2,name2,name1,name1,name2,geneCount.sum()),ha="center",wrap=True)
    plt.tight_layout()
    with open(legendOut,'w') as outFile:
        start_rtf(outFile)
        outFile.write(r'\b Figure. Number of transcripts per gene comparison\b0 \line Percent of genes with equal number of transcripts between {} and {} (Match, blue), more transcripts in {} compared to {} ({} Greater, dark orange), more transcripts in {} compared to {} ({} Greater, gray), and genes exclusive to {} ({} Only, light orange) or {} ({} Only, pink). {} genes total. Transcriptome comparisons performed by TranD [1].'.format(
            name1,name2,name1,name2,name1,name2,name1,name2,name1,name1,name2,name2,int(geneCount.sum())))
        outFile.write(r' \line \line 1. Nanni A., et al. (2021). TranD: Transcript Distance a precise nucleotide level comparison of transcript models for long reads. github')
        end_rtf(outFile)

def plot_transcript_in_gene_split_pie(md_data,f1_odds,f2_odds,name1,name2,legendOut):
    # Get counts
    geneF1only = f1_odds.groupby('gene_id')['transcript_id'].nunique().reset_index().rename(columns={'transcript_id':'num_transcript_in_gene_'+name1})
    geneF1only['num_transcript_in_gene_'+name2] = 0
    geneF1only['transcript_in_gene'] = name1+'_only'
    geneF2only = f2_odds.groupby('gene_id')['transcript_id'].nunique().reset_index().rename(columns={'transcript_id':'num_transcript_in_gene_'+name2})
    geneF2only['num_transcript_in_gene_'+name1] = 0
    geneF2only['transcript_in_gene'] = name2+'_only'

    # Make dictionary of names for plot and order
    geneCountNameDict = {'match':'Match',name1+'_greater':name1+" Greater",name2+'_greater':name2+" Greater",name1+'_only':name1+" Only",name2+'_only':name2+" Only"}
    geneCountOrderDict = {'match':0,name1+'_greater':1,name2+'_greater':2,name1+'_only':3,name2+'_only':4}

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
            totalGenes = geneCount.sum()
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
#    fig.text(0.5,-0.1,"Percent of genes with equal number of transcripts between {} and {} (blue),\nmore transcripts in {} compared to {} (dark orange), more transcripts in {} compared to {} (gray),\nand genes exclusive to {} (light orange) or {} (pink). {} genes total.".format(
#            name1,name2,name1,name2,name2,name1,name1,name2,totalGenes),ha="center",wrap=True)
    plt.tight_layout()
    with open(legendOut,'w') as outFile:
        start_rtf(outFile)
        outFile.write(r'\b Figure. Number of transcripts per gene comparison\b0 \line Percent of genes with equal number of transcripts between {} and {} (Match, blue), more transcripts in {} compared to {} ({} Greater, dark orange), more transcripts in {} compared to {} ({} Greater, gray), and genes exclusive to {} ({} Only, light orange) or {} ({} Only, pink). {} genes total. Transcriptome comparisons performed by TranD [1].'.format(
            name1,name2,name1,name2,name1,name2,name1,name2,name1,name1,name2,name2,int(totalGenes)))
        outFile.write(r' \line \line 1. Nanni A., et al. (2021). TranD: Transcript Distance a precise nucleotide level comparison of transcript models for long reads. github')
        end_rtf(outFile)

def plot_gene_stack(md_data,name1,name2,legendOut,useProp=False):
    genePairDF = md_data.groupby(['transcript_in_gene','recip_min_pair_in_gene'])['gene_id'].nunique().reset_index().pivot_table(index=['transcript_in_gene'],columns='recip_min_pair_in_gene',values='gene_id')
    for col in ['reciprocal_pairs','partial_reciprocal_pairs','no_reciprocal_pairs']:
        if col not in genePairDF.columns:
            genePairDF[col] = 0
    genePairDF = genePairDF.fillna(0)
    genePairDF = genePairDF.rename(columns={'reciprocal_pairs':'Reciprocal Pairs',
                                            'partial_reciprocal_pairs':'Partial Reciprocal Pairs',
                                            'no_reciprocal_pairs':'No Reciprocal Pairs'},
                                    index={name1+'_greater':name1+" Greater",
                                           name2+'_greater':name2+" Greater",
                                           "match":"Match"})
    if useProp:
        legendText = "Within genes with equal numbers of transcript in {} and {} (Match), the proportion of genes where all transcripts have reciprocal minimum matches (Match:Reciprocal Pairs, n = {}), at least one but not all transcript pairs are reciprocal minimum matches (Match:Partial Reciprocal Pairs, n = {}), and no pairs are recirpocal minimum matches (Match:No Reciprocal Pairs, n = {}).".format(
                name1,name2,int(genePairDF.loc['Match','Reciprocal Pairs']),int(genePairDF.loc['Match','Partial Reciprocal Pairs']),int(genePairDF.loc['Match','No Reciprocal Pairs']))
        legendText = legendText + "Within genes with more transcripts in {} compared to {} ({} Greater), the proportion of genes where all transcripts in {} have reciprocal minimum matches to a subset of {} ({} Greater:Reciprocal Pairs, n = {}), at least one but not all transcript pairs are reciprocal minimum matches ({} Greater:Partial Reciprocal Pairs, n = {}), and no pairs are recirpocal minimum matches ({} Greater:No Reciprocal Pairs, n = {}).".format(
                name1,name2,name1,name2,name1,name1,int(genePairDF.loc[name1+" Greater",'Reciprocal Pairs']),name1,int(genePairDF.loc[name1+" Greater",'Partial Reciprocal Pairs']),name1,int(genePairDF.loc[name1+" Greater",'No Reciprocal Pairs']))
        legendText = legendText + "Within genes with more transcripts in {} compared to {} ({} Greater), the proportion of genes where all transcripts in {} have reciprocal minimum matches to a subset of {} ({} Greater:Reciprocal Pairs, n = {}), at least one but not all transcript pairs are reciprocal minimum matches ({} Greater:Partial Reciprocal Pairs, n = {}), and no pairs are recirpocal minimum matches ({} Greater:No Reciprocal Pairs, n = {}).".format(
                name2,name1,name2,name1,name2,name2,int(genePairDF.loc[name2+" Greater",'Reciprocal Pairs']),name2,int(genePairDF.loc[name2+" Greater",'Partial Reciprocal Pairs']),name2,int(genePairDF.loc[name2+" Greater",'No Reciprocal Pairs']))
        genePairDF = genePairDF.div(genePairDF.sum(axis=1), axis=0)
        title = "Proportion of Genes"
    else:
        legendText = "Within genes with equal numbers of transcript in {} and {} (Match), the number of genes where all transcripts have reciprocal minimum matches (Match:Reciprocal Pairs, n = {}), at least one but not all transcript pairs are reciprocal minimum matches (Match:Partial Reciprocal Pairs, n = {}), and no pairs are recirpocal minimum matches (Match:No Reciprocal Pairs, n = {}).".format(
                name1,name2,int(genePairDF.loc['Match','Reciprocal Pairs']),int(genePairDF.loc['Match','Partial Reciprocal Pairs']),int(genePairDF.loc['Match','No Reciprocal Pairs']))
        legendText = legendText + " Within genes with more transcripts in {} compared to {} ({} Greater), the number of genes where all transcripts in {} have reciprocal minimum matches to a subset of {} ({} Greater:Reciprocal Pairs, n = {}), at least one but not all transcript pairs are reciprocal minimum matches ({} Greater:Partial Reciprocal Pairs, n = {}), and no pairs are recirpocal minimum matches ({} Greater:No Reciprocal Pairs, n = {}).".format(
                name1,name2,name1,name2,name1,name1,int(genePairDF.loc[name1+" Greater",'Reciprocal Pairs']),name1,int(genePairDF.loc[name1+" Greater",'Partial Reciprocal Pairs']),name1,int(genePairDF.loc[name1+" Greater",'No Reciprocal Pairs']))
        legendText = legendText + " Within genes with more transcripts in {} compared to {} ({} Greater), the number of genes where all transcripts in {} have reciprocal minimum matches to a subset of {} ({} Greater:Reciprocal Pairs, n = {}), at least one but not all transcript pairs are reciprocal minimum matches ({} Greater:Partial Reciprocal Pairs, n = {}), and no pairs are recirpocal minimum matches ({} Greater:No Reciprocal Pairs, n = {}).".format(
                name2,name1,name2,name1,name2,name2,int(genePairDF.loc[name2+" Greater",'Reciprocal Pairs']),name2,int(genePairDF.loc[name2+" Greater",'Partial Reciprocal Pairs']),name2,int(genePairDF.loc[name2+" Greater",'No Reciprocal Pairs']))
        title = "Number of Genes"
    genePairDF.plot(kind='bar',figsize=(10,10),stacked=True,rot=45,colormap=ListedColormap(sns.color_palette('colorblind').as_hex()))
    plt.ylabel(title)
    plt.xlabel("")
    plt.legend(bbox_to_anchor=(1.04,1), loc="upper left")
#    plt.text(1,-0.15 * genePairDF.sum(axis=1).max(),legendText,ha="center",va="top",wrap=True)
    plt.tight_layout()
    with open(legendOut,'w') as outFile:
        start_rtf(outFile)
        outFile.write(r'\b Figure. Reciprocal pairs per gene \b0 \line {} Transcriptome comparisons performed by TranD [1].'.format(legendText))
        outFile.write(r' \line \line 1. Nanni A., et al. (2021). TranD: Transcript Distance a precise nucleotide level comparison of transcript models for long reads. github')
        end_rtf(outFile)


def get_pie_label(pct, allvals):
    absolute = int(round(pct/100.*np.sum(allvals)))
    return "{:.1f}%\n({:d})".format(pct, absolute)

def plot_transcript_in_gene_upset(md_data,f1_odds,f2_odds,name1,name2,legendOut):
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
    legendText = "Number of genes with the specified number of transcripts in {} and {} indicated by the black dots below the histogram of genes counts. Columns with a single black dot represent the genes exclusive to {} (n = {}) or {} (n = {}). Genes with more than one dot are in both {} and {} (n = {})".format(name1,name2,name1,geneAll['transcript_in_gene'].value_counts()[name1+'_only'],name2,geneAll['transcript_in_gene'].value_counts()[name2+'_only'],name1,name2,geneAll['transcript_in_gene'].value_counts()[['match',name1+'_greater',name2+'_greater']].sum())
    plot_upset(geneAll.set_index(colList),"Number of Genes with Each Number of Transcripts")
#    plt.text(0.5*len(geneAll.groupby(colList)['gene_id'].count()),-0.2*geneAll.groupby(colList)['gene_id'].count().max()*len(colList),legendText,ha="center",va="top",wrap=True)
    with open(legendOut,'w') as outFile:
        start_rtf(outFile)
        outFile.write(r'\b Figure. Number of transcripts per gene comparison UpSet plot \b0 \line {} Transcriptome comparisons performed by TranD [1].'.format(legendText))
        outFile.write(r' \line \line 1. Nanni A., et al. (2021). TranD: Transcript Distance a precise nucleotide level comparison of transcript models for long reads. github.')
        end_rtf(outFile)
        
def plot_recip_min_pair_AS_upset(md_data,name1,name2,legendOut):
    recipMinPairAS = md_data[md_data['flag_recip_min_match']==1][['transcript_1','transcript_2',
                            'flag_alt_exon_recip_min_match','flag_alt_donor_acceptor_recip_min_match',
                            'flag_IR_recip_min_match','flag_5_variation_recip_min_match',
                            'flag_3_variation_recip_min_match','flag_no_shared_nt_recip_min_match']].copy()
    recipMinPairAS = recipMinPairAS.rename(columns={'flag_alt_exon_recip_min_match':'Alt. Exon',
                                                    'flag_alt_donor_acceptor_recip_min_match':'Alt. Donor/Acceptor',
                                                    'flag_IR_recip_min_match':'Intron Retention',
                                                    'flag_5_variation_recip_min_match':'5\' Variation',
                                                    'flag_3_variation_recip_min_match':'3\' Variation',
                                                    'flag_no_shared_nt_recip_min_match':'No Shared Nucleotides'})
    plot_upset(recipMinPairAS.set_index([c for c in recipMinPairAS.columns if "transcript" not in c]),
               "Number of Reciprocal Minimum Pairs with AS Categories")
    legendText = "Number of reciprocal minimum pairs with the specified types of alternative splicing between {} and {} indicated by the black dots below the histogram of pair counts (n = {} total pairs). Pairs with no shared nucleotides are in the same genes but the coordinates do not overlap.".format(name1,name2,len(recipMinPairAS))
#    plt.text(0.5*len(recipMinPairAS.groupby([c for c in recipMinPairAS.columns if "transcript" not in c])['transcript_1'].count()),-0.2*recipMinPairAS.groupby([c for c in recipMinPairAS.columns if "transcript" not in c])['transcript_1'].count().max()*len([c for c in recipMinPairAS.columns if "transcript" not in c]),legendText,ha="center",va="top",wrap=True)
    with open(legendOut,'w') as outFile:
        start_rtf(outFile)
        outFile.write(r'\b Figure. Reciprocal minimum pair alternative splicing \b0 \line {} Transcriptome comparisons performed by TranD [1].'.format(legendText))
        outFile.write(r' \line \line 1. Nanni A., et al. (2021). TranD: Transcript Distance a precise nucleotide level comparison of transcript models for long reads. github.')
        end_rtf(outFile)
        
def plot_recip_min_pair_AS_upset_nt_box(md_data,name1,name2,legendOut):
    recipMinPairAS = md_data[md_data['flag_recip_min_match']==1][['transcript_1','transcript_2',
                            'flag_alt_exon_recip_min_match','flag_alt_donor_acceptor_recip_min_match',
                            'flag_IR_recip_min_match','flag_5_variation_recip_min_match',
                            'flag_3_variation_recip_min_match','flag_no_shared_nt_recip_min_match',
                            'num_nt_diff','prop_nt_diff']].copy()
    recipMinPairAS = recipMinPairAS.rename(columns={'flag_alt_exon_recip_min_match':'Alt. Exon',
                                                    'flag_alt_donor_acceptor_recip_min_match':'Alt. Donor/Acceptor',
                                                    'flag_IR_recip_min_match':'Intron Retention',
                                                    'flag_5_variation_recip_min_match':'5\' Variation',
                                                    'flag_3_variation_recip_min_match':'3\' Variation',
                                                    'flag_no_shared_nt_recip_min_match':'No Shared Nucleotides',
                                                    'num_nt_diff':'Number NT Different',
                                                    'prop_nt_diff':'Proportion NT Different'})
    plot_upset(recipMinPairAS.set_index([c for c in recipMinPairAS.columns if ("transcript" not in c) and ("NT" not in c)]),
               "Number of Reciprocal Minimum Pairs with AS Categories",['Number NT Different','Proportion NT Different'])
    legendText = "Number of reciprocal minimum pairs with the specified types of alternative splicing between {} and {} indicated by the black dots below the histogram of pair counts (n = {} total pairs). Pairs with no shared nucleotides are in the same genes but with nonoverlapping coordinates. Box plots of the number (blue) and proportion (orange) of nucleotide (NT) differences between the pairs represented in the histogram.".format(name1,name2,len(recipMinPairAS))
#    plt.text(0.5*len(recipMinPairAS.groupby([c for c in recipMinPairAS.columns if ("transcript" not in c) and ("NT" not in c)])['transcript_1'].count()),-5.25,legendText,ha="center",va="top",wrap=True)
    with open(legendOut,'w') as outFile:
        start_rtf(outFile)
        outFile.write(r'\b Figure. Reciprocal minimum pair alternative splicing \b0 \line {} Transcriptome comparisons performed by TranD [1].'.format(legendText))
        outFile.write(r' \line \line 1. Nanni A., et al. (2021). TranD: Transcript Distance a precise nucleotide level comparison of transcript models for long reads. github.')
        end_rtf(outFile)
        
def plot_gene_AS_upset(md_data,name1,name2,legendOut):
    geneRecipMatchAS = md_data.groupby('gene_id')[['flag_alt_exon_recip_min_match','flag_alt_donor_acceptor_recip_min_match',
                                                   'flag_IR_recip_min_match','flag_5_variation_recip_min_match',
                                                   'flag_3_variation_recip_min_match',
                                                   'flag_no_shared_nt_recip_min_match']].max().astype(bool).reset_index()
    geneRecipMatchAS = geneRecipMatchAS.rename(columns={'flag_alt_exon_recip_min_match':'Alt. Exon',
                                                        'flag_alt_donor_acceptor_recip_min_match':'Alt. Donor/Acceptor',
                                                        'flag_IR_recip_min_match':'Intron Retention',
                                                        'flag_5_variation_recip_min_match':'5\' Variation',
                                                        'flag_3_variation_recip_min_match':'3\' Variation',
                                                        'flag_no_shared_nt_recip_min_match':'No Shared Nucleotides'})
    plot_upset(geneRecipMatchAS.set_index([c for c in geneRecipMatchAS.columns if c != "gene_id"]),
               "Number of Genes with AS Categories in Reciprocal Minimum Pairs")
    legendText = "Number of genes with the specified types of alternative splicing in only reciprocal minimum pairs between {} and {} indicated by the black dots below the histogram of gene counts (n = {} genes with {} reciprocal minimum pairs). Genes with \"No Shared Nucleotides\" have a pair of transcripts with nonoverlapping coordinates.".format(name1,name2,len(geneRecipMatchAS),len(md_data[md_data['flag_recip_min_match']==1]))
#    plt.text(0.5*len(geneRecipMatchAS.groupby([c for c in geneRecipMatchAS.columns if c != "gene_id"])['gene_id'].count()),-0.2*geneRecipMatchAS.groupby([c for c in geneRecipMatchAS.columns if c != "gene_id"])['gene_id'].count().max()*len([c for c in geneRecipMatchAS.columns if c != "gene_id"]),legendText,ha="center",va="top",wrap=True)
    with open(legendOut,'w') as outFile:
        start_rtf(outFile)
        outFile.write(r'\b Figure. Alternative splicing in reciprocal minimum pairs of genes \b0 \line {} Transcriptome comparisons performed by TranD [1].'.format(legendText))
        outFile.write(r' \line \line 1. Nanni A., et al. (2021). TranD: Transcript Distance a precise nucleotide level comparison of transcript models for long reads. github.')
        end_rtf(outFile)
        
def plot_gene_AS_upset_nt_box(md_data,name1,name2,legendOut):
    recipMinPairAS = md_data[md_data['flag_recip_min_match']==1][['gene_id',
                            'flag_alt_exon_recip_min_match','flag_alt_donor_acceptor_recip_min_match',
                            'flag_IR_recip_min_match','flag_5_variation_recip_min_match',
                            'flag_3_variation_recip_min_match','flag_no_shared_nt_recip_min_match',
                            'num_recip_min_match_in_gene','num_transcript_in_gene_'+name1,
                            'num_transcript_in_gene_'+name2,'num_nt_diff','prop_nt_diff']].copy()
    # Ensure number of nt different are int and float values
    recipMinPairAS['num_nt_diff'] = recipMinPairAS['num_nt_diff'].astype(int)
    recipMinPairAS['prop_nt_diff'] = recipMinPairAS['num_nt_diff'].astype(float)
    geneRecipMatchAS = recipMinPairAS.groupby('gene_id').agg({'flag_alt_exon_recip_min_match':'max',
                                                              'flag_alt_donor_acceptor_recip_min_match':'max',
                                                              'flag_IR_recip_min_match':'max',
                                                              'flag_5_variation_recip_min_match':'max',
                                                              'flag_3_variation_recip_min_match':'max',
                                                              'flag_no_shared_nt_recip_min_match':'max',
                                                              'num_recip_min_match_in_gene':'max',
                                                              'num_transcript_in_gene_'+name1:'max',
                                                              'num_transcript_in_gene_'+name2:'max',
                                                              'num_nt_diff':'mean',
                                                              'prop_nt_diff':'mean'}).reset_index()
    geneFlagCols = ['flag_alt_exon_recip_min_match','flag_alt_donor_acceptor_recip_min_match',
                    'flag_IR_recip_min_match','flag_5_variation_recip_min_match',
                    'flag_3_variation_recip_min_match','flag_no_shared_nt_recip_min_match']
    geneRecipMatchAS[geneFlagCols] = geneRecipMatchAS[geneFlagCols].astype(bool)
    geneRecipMatchAS = geneRecipMatchAS.rename(columns={'flag_alt_exon_recip_min_match':'Alt. Exon',
                                                        'flag_alt_donor_acceptor_recip_min_match':'Alt. Donor/Acceptor',
                                                        'flag_IR_recip_min_match':'Intron Retention',
                                                        'flag_5_variation_recip_min_match':'5\' Variation',
                                                        'flag_3_variation_recip_min_match':'3\' Variation',
                                                        'flag_no_shared_nt_recip_min_match':'No Shared Nucleotides',
                                                        'num_recip_min_match_in_gene':'# Recip. Min.\nMatch Transcripts',
                                                        'num_transcript_in_gene_'+name1:'# Transcripts\nin '+name1,
                                                        'num_transcript_in_gene_'+name2:'# Transcripts\nin '+name2,
                                                        'num_nt_diff':'Avg #\nNT Different',
                                                        'prop_nt_diff':'Avg Proportion\nNT Different'})
    plot_upset(geneRecipMatchAS.set_index([c for c in geneRecipMatchAS.columns if ("gene" not in c) and ("NT" not in c) and ("Transcripts" not in c) ]),
               "",['# Recip. Min.\nMatch Transcripts','# Transcripts\nin '+name1,'# Transcripts\nin '+name2,'Avg #\nNT Different','Avg Proportion\nNT Different'])
    legendText = "Number of genes with the specified types of alternative splicing in only reciprocal minimum pairs between {} and {} indicated by the black dots below the histogram of gene counts (n = {} genes with {} reciprocal minimum pairs). Box plots represent the number of reciprocal minimum pairs (blue), number of transcripts in {} (orange) and {} (green), and the average number (brown) and proportion (purple) of nucleotides different between the pairs. Genes with \"No Shared Nucleotides\" have a pair of transcripts with nonoverlapping coordinates.".format(name1,name2,len(geneRecipMatchAS),len(md_data[md_data['flag_recip_min_match']==1]),name1,name2)
#    plt.text(0.5*len(geneRecipMatchAS.groupby([c for c in geneRecipMatchAS.columns if ("gene" not in c) and ("NT" not in c) and ("Transcripts" not in c) ])['gene_id'].count()),-8.25*geneRecipMatchAS.groupby([c for c in geneRecipMatchAS.columns if ("gene" not in c) and ("NT" not in c) and ("Transcripts" not in c) ])['Avg Proportion\nNT Different'].max().max(),legendText,ha="center",va="top",wrap=True)
    with open(legendOut,'w') as outFile:
        start_rtf(outFile)
        outFile.write(r'\b Figure. Alternative splicing in reciprocal minimum pairs of genes \b0 \line {} Transcriptome comparisons performed by TranD [1].'.format(legendText))
        outFile.write(r' \line \line 1. Nanni A., et al. (2021). TranD: Transcript Distance a precise nucleotide level comparison of transcript models for long reads. github.')
        end_rtf(outFile)
        
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

def plot_gene_avg_nt_diff_pairs(md_data,name1,name2,legendOut,zoomMean=False):
    minPairNT = md_data[['gene_id','transcript_in_gene','recip_min_pair_in_gene',
                         'flag_recip_min_match','num_nt_diff']].copy()
    # Ensure number of nt different are int and float values
    minPairNT['num_nt_diff'] = minPairNT['num_nt_diff'].astype(int)
    minPairNT['prop_nt_diff'] = minPairNT['num_nt_diff'].astype(float)
    minPairNT['num_nt_diff_recip_min'] = np.where(minPairNT['flag_recip_min_match']==1,minPairNT['num_nt_diff'],0)
    minPairNT['num_nt_diff_extra'] = np.where(minPairNT['flag_recip_min_match']==0,minPairNT['num_nt_diff'],0)
    minPairNTGene = minPairNT.groupby('gene_id').agg({'transcript_in_gene':'first',
                                                'recip_min_pair_in_gene':'first',
                                                'num_nt_diff_recip_min':'mean',
                                                'num_nt_diff_extra':'mean'}).reset_index()
    recipConditions = [(minPairNTGene['transcript_in_gene']=="match")&(minPairNTGene['recip_min_pair_in_gene']=="reciprocal_pairs"),
                       (minPairNTGene['transcript_in_gene']=="match")&(minPairNTGene['recip_min_pair_in_gene']=="partial_reciprocal_pairs"),
                       (minPairNTGene['transcript_in_gene']=="match")&(minPairNTGene['recip_min_pair_in_gene']=="no_reciprocal_pairs"),
                       (minPairNTGene['transcript_in_gene']==name1+'_greater')&(minPairNTGene['recip_min_pair_in_gene']=="reciprocal_pairs"),
                       (minPairNTGene['transcript_in_gene']==name1+'_greater')&(minPairNTGene['recip_min_pair_in_gene']=="partial_reciprocal_pairs"),
                       (minPairNTGene['transcript_in_gene']==name1+'_greater')&(minPairNTGene['recip_min_pair_in_gene']=="no_reciprocal_pairs"),
                       (minPairNTGene['transcript_in_gene']==name2+'_greater')&(minPairNTGene['recip_min_pair_in_gene']=="reciprocal_pairs"),
                       (minPairNTGene['transcript_in_gene']==name2+'_greater')&(minPairNTGene['recip_min_pair_in_gene']=="partial_reciprocal_pairs"),
                       (minPairNTGene['transcript_in_gene']==name2+'_greater')&(minPairNTGene['recip_min_pair_in_gene']=="no_reciprocal_pairs")]
    recipChoices = ["Match:Reciprocal Pairs","Match:Partial Reciprocal Pairs","Match:No Reciprocal Pairs",
                    name1+" Greater:Reciprocal Pairs",name1+" Greater:Partial Reciprocal Pairs",name1+" Greater:No Reciprocal Pairs",
                    name2+" Greater:Reciprocal Pairs",name2+" Greater:Partial Reciprocal Pairs",name2+" Greater:No Reciprocal Pairs"]
    minPairNTGene['Category'] = np.select(recipConditions,recipChoices,"NC")
    minPairGeneOrder = {"Match:Reciprocal Pairs":0,"Match:Partial Reciprocal Pairs":1,"Match:No Reciprocal Pairs":2,
                        name1+" Greater:Reciprocal Pairs":3,name1+" Greater:Partial Reciprocal Pairs":4,name1+" Greater:No Reciprocal Pairs":5,
                        name2+" Greater:Reciprocal Pairs":6,name2+" Greater:Partial Reciprocal Pairs":7,name2+" Greater:No Reciprocal Pairs":8}
    minPairNTGene['Category'] = pd.Categorical(minPairNTGene['Category'],categories=sorted(minPairGeneOrder,key=minPairGeneOrder.get),ordered=True)
    minPairNTGene = minPairNTGene.sort_values('Category')
    colorPalleteRecip = {"Match:Reciprocal Pairs":"#0173b2","Match:Partial Reciprocal Pairs":"#56b4e9","Match:No Reciprocal Pairs":"#029e73",
                    name1+" Greater:Reciprocal Pairs":"#d55e00",name1+" Greater:Partial Reciprocal Pairs":"#de8f05",name1+" Greater:No Reciprocal Pairs":"#ece133",
                    name2+" Greater:Reciprocal Pairs":"#949494",name2+" Greater:Partial Reciprocal Pairs":"#fbafe4",name2+" Greater:No Reciprocal Pairs":"#ca9161"}
    if not zoomMean:
        limit = minPairNTGene[['num_nt_diff_recip_min','num_nt_diff_extra']].max().max()
    else:
        limit=round(max(minPairNTGene[minPairNTGene['num_nt_diff_recip_min']>0]['num_nt_diff_recip_min'].mean(),
                        minPairNTGene[minPairNTGene['num_nt_diff_extra']>0]['num_nt_diff_extra'].mean()))
    plt.figure(figsize=(8,5))
    sns.scatterplot(data=minPairNTGene,x='num_nt_diff_recip_min',
                    y='num_nt_diff_extra',hue='Category',palette=colorPalleteRecip)
    plt.xlim(0,limit)
    plt.ylim(0,limit)
    plt.xlabel("Avg. # NT Different in Recip. Min. Pairs")
    plt.ylabel("Avg. # NT Different in Min. Pairs of Extras")
    plt.legend(bbox_to_anchor=(1.04,1), loc="upper left")
    plt.plot((0,limit),(0,limit),color='r')
    legendText = "Each point represents a gene plotted by the average number of nucleotides different between reciprocal minimum pairs (X-axis) and minimum pairs of extras without a reciprocal minimum pair (Y-axis). Genes are colored by categorizations of the comparison between {} and {}".format(name1,name2)
#    plt.text(0.5*limit,-0.15*limit,legendText,ha="center",va="top",wrap=True)
    plt.tight_layout()
    with open(legendOut,'w') as outFile:
        start_rtf(outFile)
        outFile.write(r'\b Figure. Number of nucleotides different between reciprocal and nonreciprocal minimum pairs \b0 \line {} Transcriptome comparisons performed by TranD [1].'.format(legendText))
        outFile.write(r' \line \line 1. Nanni A., et al. (2021). TranD: Transcript Distance a precise nucleotide level comparison of transcript models for long reads. github.')
        end_rtf(outFile)
        
def start_rtf(outFile):
    outFile.write(r'{\rtf1\ansi\ansicpg1252\deff0\deflang1033{\fonttbl{\f0\fswiss\fcharset0 Arial;}}')

def end_rtf(outFile):
    outFile.write(" \line \line \line \i Disclaimer:  While automated captions of TranD have been carefully constructed, users are advised to verify caption contents before use and report any errors to the TranD github.\i0 ")
    outFile.write(r'}\n\x00')