#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 13 15:03:09 2021

@author: adalena.nanni
"""

import pandas as pd
import numpy as np

td_df_cols = ['gene_id','transcript_1','transcript_2','num_junction_T1_only','num_junction_T2_only',
             'num_junction_shared','prop_junction_diff','prop_junction_similar',
             'junction_T1_only','junction_T2_only','junction_shared','num_ER_T1_only',
             'num_ER_T2_only','num_ER_shared','prop_ER_diff','prop_ER_similar',
             'ER_T1_only','ER_T2_only','ER_shared','num_fragment_T1_only',
             'num_fragment_T2_only','num_fragment_shared','prop_fragment_diff',
             'prop_fragment_similar','fragment_T1_only','fragment_T2_only','fragment_shared',
             'num_fragment_singletons_T1_only','num_fragment_singletons_T2_only',
             'num_fragment_singletons_shared','num_IR_fragment_T1','num_IR_fragment_T2',
             'IR_fragment_T1','IR_fragment_T2','num_nt_shared','num_nt_T1_only',
             'num_nt_T2_only','total_nt','prop_nt_diff','prop_nt_similar','num_nt_T1_only_in_shared_ER',
             'num_nt_T2_only_in_shared_ER','num_nt_shared_in_shared_ER','total_nt_in_shared_ER',
             'prop_nt_diff_in_shared_ER','prop_nt_similar_in_shared_ER','num_nt_T1_only_in_unique_ER',
             'num_nt_T2_only_in_unique_ER']


def calculate_distance(out_df,junction_df,gene_id,tx1_name,tx2_name,fsm=False):
    # Make empty Series to place distance values
    singlePair = pd.Series(index=td_df_cols,dtype=object)
    singlePair[['gene_id','transcript_1','transcript_2']] = gene_id,tx1_name,tx2_name

    # Set exon fragment and exon region ID values based on genomic coordinates
    out_df['fragment_id'] = out_df['ef_chr'].map(str)+":"+out_df['ef_start'].map(str)+":"+out_df['ef_end'].map(str)+":"+out_df['ef_strand'].map(str)
    out_df['region_id'] = out_df['er_chr'].map(str)+":"+out_df['er_start'].map(str)+":"+out_df['er_end'].map(str)+":"+out_df['er_strand'].map(str)

    # Flag shared fragments
    out_df['flag_shared'] = np.where((out_df['transcript_id']==tx1_name+"|"+tx2_name)|(out_df['transcript_id']==tx2_name+"|"+tx1_name),1,0)

    # Flag singleton fragments (where the fragment is a full exon region)
    out_df['num_fragment_in_ER'] = out_df.groupby('er_id')['ef_id'].transform('count')
    out_df['flag_singleton'] = np.where(out_df['num_fragment_in_ER']==1,1,0)
    
    # Add fragment lengths (end - start)
    out_df['fragment_length'] = out_df['ef_end'].map(int) - out_df['ef_start'].map(int)

    # Get distance measures for junctions, exon regions (ER),and exon fragments (EF)
    singlePair = get_junction_distance(singlePair,junction_df,tx1_name,tx2_name,fsm=fsm)
    singlePair, ERSharedSet = get_ER_distance(singlePair,out_df,tx1_name,tx2_name,fsm=fsm)
    singlePair = get_EF_distance(singlePair,out_df,tx1_name,tx2_name,ERSharedSet)
    
    # Set flags for different alternative splicing (AS) events
    singlePair = set_AS_flags(singlePair)
    
    # Return distance of transcript pair
    return singlePair


def get_junction_distance(singlePair,junction_df,tx1_name,tx2_name,fsm=False):
    # Check if transcript pair shares all junctions (FSM, full-splice match)
    if fsm:
        # Check if both transcripts are monoexon (no junctions)
        if len(junction_df) == 0:
            singlePair[['num_junction_T1_only','num_junction_T2_only','num_junction_shared']] = 0
            singlePair['prop_junction_diff'] = 0
            singlePair['prop_junction_similar'] = 1
            singlePair[['junction_T1_only','junction_T2_only','junction_shared']] = ""
        # FSM transcripts are multiexon
        else:
            singlePair[['num_junction_T1_only','num_junction_T2_only']] = 0
            singlePair['num_junction_shared'] = len(junction_df[junction_df['transcript_id']==tx1_name])
            singlePair['prop_junction_diff'] = 0
            singlePair['prop_junction_similar'] = 1
            singlePair[['junction_T1_only','junction_T2_only']] = ""
            singlePair['junction_shared'] = "|".join(junction_df[junction_df['transcript_id']==tx1_name]['coords'])
    # Transcript pair does not share all junctions
    else:
        # Check if both transcripts are monoexon but do not overlap (not fsm)
        if len(junction_df) == 0:
            singlePair[['num_junction_T1_only','num_junction_T2_only','num_junction_shared']] = 0
            singlePair['prop_junction_diff'] = 0
            singlePair['prop_junction_similar'] = 1
            singlePair[['junction_T1_only','junction_T2_only','junction_shared']] = ""            
        # Check if one of the transcripts is monoexon (only junctions from one transcript and not the other)
        elif junction_df['transcript_id'].nunique() == 1:
        # Only T1 is monoexon
            if tx1_name in junction_df['transcript_id'].unique():
                singlePair[['num_junction_T1_only','num_junction_shared']] = 0
                singlePair['prop_junction_similar'] = 0
                singlePair['prop_junction_diff'] = 1
                singlePair[['junction_T1_only','junction_shared']] = ""
                singlePair['num_junction_T2_only'] = len(junction_df)
                singlePair['junction_T2_only'] = "|".join(junction_df['coords'])
        # Only T2 is monoexon
            else:
                singlePair[['num_junction_T2_only','num_junction_shared']] = 0
                singlePair['prop_junction_similar'] = 0
                singlePair['prop_junction_diff'] = 1
                singlePair[['junction_T2_only','junction_shared']] = ""
                singlePair['num_junction_T1_only'] = len(junction_df)
                singlePair['junction_T1_only'] = "|".join(junction_df['coords'])
        # Both transcripts multi-exon
        else:
            t1Junc = junction_df[junction_df['transcript_id']==tx1_name]['coords']
            t2Junc = junction_df[junction_df['transcript_id']==tx2_name]['coords']
            juncSharedSet = set(t1Junc).intersection(set(t2Junc))
            juncT1Set = set(t1Junc).difference(set(t2Junc))
            juncT2Set = set(t2Junc).difference(set(t1Junc))
            singlePair['num_junction_T1_only'] = len(juncT1Set)
            singlePair['num_junction_T2_only'] = len(juncT2Set)
            singlePair['num_junction_shared'] = len(juncSharedSet)
            singlePair['prop_junction_diff'] = (singlePair['num_junction_T1_only'] + singlePair['num_junction_T2_only'])/(singlePair['num_junction_T1_only'] + singlePair['num_junction_T2_only'] + singlePair['num_junction_shared'])
            singlePair['prop_junction_similar'] = 1 - singlePair['prop_junction_diff']
            singlePair['junction_T1_only'] = "|".join(sorted(juncT1Set, key=lambda x: t1Junc[t1Junc == x].index[0]))
            singlePair['junction_T2_only'] = "|".join(sorted(juncT2Set, key=lambda x: t2Junc[t2Junc == x].index[0]))
            singlePair['junction_shared'] = "|".join(sorted(juncSharedSet, key=lambda x: t2Junc[t2Junc == x].index[0]))
    return singlePair


def get_ER_distance(singlePair,out_df,tx1_name,tx2_name,fsm=False):
    # Check if transcript pair shares all junctions (FSM, full-splice match)
    if fsm:
        # All junctions are shared, therefore all ER are shared
        singlePair[['num_ER_T1_only','num_ER_T2_only']] = 0
        singlePair['num_ER_shared'] = out_df['er_id'].nunique()
        singlePair['prop_ER_diff'] = 0
        singlePair['prop_ER_similar'] = 1
        singlePair[['ER_T1_only','ER_T2_only']] = ""
        singlePair['ER_shared'] = "|".join(out_df['region_id'].unique())
        ERSharedSet = out_df['region_id'].drop_duplicates()
    else:
        # Check for shared and unique ER
        ERall = out_df.groupby('region_id').agg({'transcript_id':(lambda x: max(x,key=len)),'flag_shared':'max'}).reset_index()
        ERSharedSet = ERall[ERall['flag_shared']==1]['region_id']
        ERT1Set = ERall[(ERall['flag_shared']==0)&(ERall['transcript_id']==tx1_name)]['region_id']
        ERT2Set = ERall[(ERall['flag_shared']==0)&(ERall['transcript_id']==tx2_name)]['region_id']
        singlePair['num_ER_T1_only'] = len(ERT1Set)
        singlePair['num_ER_T2_only'] = len(ERT2Set)
        singlePair['num_ER_shared'] = len(ERSharedSet)
        singlePair['prop_ER_diff'] = (singlePair['num_ER_T1_only'] + singlePair['num_ER_T2_only'])/(singlePair['num_ER_T1_only'] + singlePair['num_ER_T2_only'] + singlePair['num_ER_shared'])
        singlePair['prop_ER_similar'] = 1 - singlePair['prop_ER_diff']
        singlePair['ER_T1_only'] = "|".join(ERT1Set)
        singlePair['ER_T2_only'] = "|".join(ERT2Set)
        singlePair['ER_shared'] = "|".join(ERSharedSet)
    return singlePair, ERSharedSet


def get_EF_distance(singlePair,out_df,tx1_name,tx2_name,ERSharedSet):
    num_shared = out_df['flag_shared'].sum()
    # Check if all fragments are shared (transcripts are identical)
    if num_shared == len(out_df):
        # All EF shared, transcripts are identical
        singlePair[['num_fragment_T1_only','num_fragment_T2_only']] = 0
        singlePair['num_fragment_shared'] = len(out_df)
        singlePair['prop_fragment_diff'] = 0
        singlePair['prop_fragment_similar'] = 1
        singlePair[['fragment_T1_only','fragment_T2_only']] = ""
        singlePair['fragment_shared'] = "|".join(out_df['fragment_id'])
        singlePair[['num_fragment_singletons_T1_only','num_fragment_singletons_T2_only']] = 0
        singlePair['num_fragment_singletons_shared'] = singlePair['num_fragment_shared']
        
        # Nucleotide differences are 0, all are shared and are in shared ER
        singlePair[['total_nt','total_nt_in_shared_ER','num_nt_shared','num_nt_shared_in_shared_ER']] = out_df['fragment_length'].sum()
        singlePair[['num_nt_T1_only','num_nt_T1_only_in_shared_ER','num_nt_T1_only_in_unique_ER',
                    'num_nt_T2_only','num_nt_T2_only_in_shared_ER','num_nt_T2_only_in_unique_ER']] = 0
        singlePair[['prop_nt_diff','prop_nt_diff_in_shared_ER']] = 0
        singlePair[['prop_nt_similar','prop_nt_similar_in_shared_ER']] = 1
        
        # No IR present between identical transcripts
        singlePair[['num_IR_fragment_T1','num_IR_fragment_T2']] = 0
        singlePair[['IR_fragment_T1','IR_fragment_T2']] = ""   
    else:
        # Not all EF shared
        fragSharedSet = out_df[out_df['flag_shared']==1]
        fragSharedSingSet = out_df[(out_df['flag_shared']==1)&(out_df['flag_singleton']==1)]
        fragT1Set = out_df[(out_df['flag_shared']==0)&(out_df['transcript_id']==tx1_name)]
        fragT1SingSet = out_df[(out_df['flag_shared']==0)&(out_df['transcript_id']==tx1_name)&(out_df['flag_singleton']==1)]
        fragT2Set = out_df[(out_df['flag_shared']==0)&(out_df['transcript_id']==tx2_name)]
        fragT2SingSet = out_df[(out_df['flag_shared']==0)&(out_df['transcript_id']==tx2_name)&(out_df['flag_singleton']==1)]
        singlePair['num_fragment_T1_only'] = len(fragT1Set)
        singlePair['num_fragment_T2_only'] = len(fragT2Set)
        singlePair['num_fragment_shared'] = len(fragSharedSet)
        singlePair['prop_fragment_diff'] = (singlePair['num_fragment_T1_only'] + singlePair['num_fragment_T2_only'])/(singlePair['num_fragment_T1_only'] + singlePair['num_fragment_T2_only'] + singlePair['num_fragment_shared'])
        singlePair['prop_fragment_similar'] = 1 - singlePair['prop_fragment_diff']
        singlePair['fragment_T1_only'] = "|".join(fragT1Set['fragment_id'])
        singlePair['fragment_T2_only'] = "|".join(fragT2Set['fragment_id'])
        singlePair['fragment_shared'] = "|".join(fragSharedSet['fragment_id'])
        singlePair['num_fragment_singletons_T1_only'] = len(fragT1SingSet)
        singlePair['num_fragment_singletons_T2_only'] = len(fragT2SingSet)
        singlePair['num_fragment_singletons_shared'] = len(fragSharedSingSet)

        # Count number of nt shared/different in all EF
        singlePair['num_nt_shared'] = fragSharedSet['fragment_length'].sum()
        singlePair['num_nt_T1_only'] = fragT1Set['fragment_length'].sum()
        singlePair['num_nt_T2_only'] = fragT2Set['fragment_length'].sum()
        singlePair['total_nt'] = singlePair['num_nt_shared'] + singlePair['num_nt_T1_only'] + singlePair['num_nt_T2_only']
        singlePair['prop_nt_diff'] = (singlePair['num_nt_T1_only'] + singlePair['num_nt_T2_only'])/(singlePair['total_nt'])
        singlePair['prop_nt_similar'] = 1 - singlePair['prop_nt_diff']
        
        # Count number of nt shared/different in EF only in shared ER
        singlePair['num_nt_T1_only_in_shared_ER'] = fragT1Set[fragT1Set['region_id'].isin(ERSharedSet)]['fragment_length'].sum()
        singlePair['num_nt_T2_only_in_shared_ER'] = fragT2Set[fragT2Set['region_id'].isin(ERSharedSet)]['fragment_length'].sum()
        singlePair['num_nt_shared_in_shared_ER'] = fragSharedSet[fragSharedSet['region_id'].isin(ERSharedSet)]['fragment_length'].sum()
        singlePair['total_nt_in_shared_ER'] = singlePair['num_nt_shared_in_shared_ER'] + singlePair['num_nt_T1_only_in_shared_ER'] + singlePair['num_nt_T2_only_in_shared_ER']
        if singlePair['total_nt_in_shared_ER'] != 0:
            singlePair['prop_nt_diff_in_shared_ER'] = (singlePair['num_nt_T1_only_in_shared_ER'] + singlePair['num_nt_T2_only_in_shared_ER'])/(singlePair['total_nt_in_shared_ER'])
            singlePair['prop_nt_similar_in_shared_ER'] = 1 - singlePair['prop_nt_diff_in_shared_ER']
        else:
            singlePair['prop_nt_diff_in_shared_ER'] = 0
            singlePair['prop_nt_similar_in_shared_ER'] = 0
        singlePair['num_nt_T1_only_in_unique_ER'] = fragT1Set[~fragT1Set['region_id'].isin(ERSharedSet)]['fragment_length'].sum()
        singlePair['num_nt_T2_only_in_unique_ER'] = fragT2Set[~fragT2Set['region_id'].isin(ERSharedSet)]['fragment_length'].sum()
        
        # Get IR distance values
        singlePair['num_IR_fragment_T1'] = len(fragT1Set[fragT1Set['ef_ir_flag'].map(int)==1])
        singlePair['num_IR_fragment_T2'] = len(fragT2Set[fragT2Set['ef_ir_flag'].map(int)==1])
        singlePair['IR_fragment_T1'] = "|".join(fragT1Set[fragT1Set['ef_ir_flag'].map(int)==1]['fragment_id'])
        singlePair['IR_fragment_T2'] = "|".join(fragT2Set[fragT2Set['ef_ir_flag'].map(int)==1]['fragment_id'])
    return singlePair


def set_AS_flags(singlePair):
    # Flag alternative splicing events
    # Alternate exons are when not all exon regions are shared
    singlePair['flag_alt_exon'] = np.where(singlePair['prop_ER_diff']>0,1,0)
    
    # Alternate 5'/3' ends are when the first/last exon regions are shared
    #   and the difference to the TSS/TTS is > 0
    # !!! Currently these are flagged as alt donor/acceptors
    # !!! Add once num_nt_diff_TSS (number of nt different between the TSS of T1 and TSS of T2, always positive or 0 if the same)
    # !!! Add once num_nt_diff_TTS (number of nt different between the TTS of T1 and TTS of T2, always positive or 0 if the same)
    # !!! Add once flag_5prime_ER_shared (1 if the 5' most ER is shared betwen T1 and T2, else 0 - using the first ER if on + strand or last ER if on - strand)
    # !!! Add once flag_3prime_ER_shared (1 if the 3' most ER is shared betwen T1 and T2, else 0 - using the last ER if on + strand or first ER if on - strand)

    # Alternate donor/acceptors in shared ER
    singlePair['flag_alt_donor_acceptor'] = np.where(singlePair['prop_nt_diff_in_shared_ER']>0,1,0)
    # !!! Will need to fix so that this is not including transcripts that share all junctions once the 5'/3' end flag set up
    # singlePair['flag_alt_donor_acceptor'] = np.where((singlePair['prop_nt_diff_in_shared_ER']>0)&(singlePair['prop_junction_diff']>0),1,0)
    # !!! May also need to fix so that nt differences in shared ER that are IR are not counted somehow

    # Intron retention 
    singlePair['flag_IR'] = np.where(singlePair['num_IR_fragment_T1']+singlePair['num_IR_fragment_T2']>0,1,0)
    
    return singlePair