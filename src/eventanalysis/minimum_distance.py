#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 20 12:49:29 2021

@author: adalena.nanni
"""

import numpy as np
from loguru import logger

md_df_cols = ['gene_id','transcript_1','transcript_2','num_junction_T1_only','num_junction_T2_only',
             'num_junction_shared','prop_junction_diff','prop_junction_similar',
             'junction_T1_only','junction_T2_only','junction_shared','num_ER_T1_only',
             'num_ER_T2_only','num_ER_shared','prop_ER_diff','prop_ER_similar',
             'ER_T1_only','ER_T2_only','ER_shared','num_fragment_T1_only',
             'num_fragment_T2_only','num_fragment_shared','prop_fragment_diff',
             'prop_fragment_similar','fragment_T1_only','fragment_T2_only','fragment_shared',
             'num_fragment_singletons_T1_only','num_fragment_singletons_T2_only',
             'num_fragment_singletons_shared','num_IR_fragment_T1','num_IR_fragment_T2',
             'IR_fragment_T1','IR_fragment_T2','num_nt_shared','num_nt_T1_only',
             'num_nt_T2_only','num_nt_diff','total_nt','prop_nt_diff','prop_nt_similar','num_nt_T1_only_in_shared_ER',
             'num_nt_T2_only_in_shared_ER','num_nt_shared_in_shared_ER','total_nt_in_shared_ER',
             'prop_nt_diff_in_shared_ER','prop_nt_similar_in_shared_ER','num_nt_T1_only_in_unique_ER',
             'num_nt_T2_only_in_unique_ER','flag_FSM','flag_IR','flag_5_variation',
             'flag_3_variation','flag_alt_donor_acceptor','flag_alt_exon','flag_nonoverlapping','num_transcript_in_gene_d1',
             'num_transcript_in_gene_d2','flag_d1_greater','flag_d2_greater','flag_equal',
             'transcript_in_gene','min_match_d1','flag_min_match_d1','min_match_d2',
             'flag_min_match_d2','flag_recip_min_match','num_recip_min_match_in_gene',
             'flag_d1_tie','flag_d2_tie','flag_identical_recip_min_match','flag_FSM_recip_min_match',
             'flag_ERM_recip_min_match','flag_ERM_noIR_recip_min_match','flag_ERM_withIR_recip_min_match',
             'num_identical_recip_min_match_in_gene','num_FSM_recip_min_match_in_gene',
             'num_ERM_recip_min_match_in_gene','num_ERM_noIR_recip_min_match_in_gene',
             'num_ERM_withIR_recip_min_match_in_gene','prop_identical_recip_min_match_in_gene',
             'prop_FSM_recip_min_match_in_gene','prop_ERM_recip_min_match_in_gene',
             'prop_ERM_noIR_recip_min_match_in_gene','prop_ERM_withIR_recip_min_match_in_gene',
             'flag_nonoverlapping_recip_min_match','flag_alt_exon_recip_min_match','flag_alt_donor_acceptor_recip_min_match',
             'flag_IR_recip_min_match','flag_5_variation_recip_min_match','flag_3_variation_recip_min_match',
             'recip_min_pair_in_gene']

def identify_min_pair(td_data, all_pairs):
    # Set names of input data
    name1 = "d1"
    name2 = "d2"
        
    # Get number of transcripts per gene
    td_data['num_transcript_in_gene_'+name1] = td_data.groupby('gene_id')['transcript_1'].transform('nunique')
    td_data['num_transcript_in_gene_'+name2] = td_data.groupby('gene_id')['transcript_2'].transform('nunique')

     # Flag groups based on transcripts per gene
    td_data['flag_'+name1+'_greater'] = np.where((td_data['num_transcript_in_gene_'+name1]>td_data['num_transcript_in_gene_'+name2])&
                                                 (td_data['num_transcript_in_gene_'+name1]>0)&(td_data['num_transcript_in_gene_'+name2]>0),1,0)
    td_data['flag_'+name2+'_greater'] = np.where((td_data['num_transcript_in_gene_'+name1]<td_data['num_transcript_in_gene_'+name2])&
                                                 (td_data['num_transcript_in_gene_'+name1]>0)&(td_data['num_transcript_in_gene_'+name2]>0),1,0)
    td_data['flag_equal'] = np.where((td_data['num_transcript_in_gene_'+name1]==td_data['num_transcript_in_gene_'+name2]),1,0)
    conditionsGene = [td_data['flag_'+name1+'_greater']==1,td_data['flag_'+name2+'_greater']==1,td_data['flag_equal']==1]
    choicesGene = [name1+"_greater",name2+"_greater","equal"]
    td_data['transcript_in_gene'] = np.select(conditionsGene,choicesGene,'missing')
    if len(td_data[td_data['transcript_in_gene']=="missing"]) > 0:
        logger.error("Unexpected variable assignment for gene transcript count comparison")
    
    # Flag the minimum pair for each transcript
    # First sort (ascending) by proportion of junctions different,
    #   proportion of ER different, and proportion of nt different
    td_data = td_data.sort_values(['prop_junction_diff','prop_ER_diff','prop_nt_diff','num_nt_diff'], ascending=True)
    td_data['min_match_'+name1] = td_data.groupby('transcript_1')['transcript_2'].transform('first')
    td_data['flag_min_match_'+name1] = np.where(td_data['min_match_'+name1]==td_data['transcript_2'],1,0)
    td_data['min_match_'+name2] = td_data.groupby('transcript_2')['transcript_1'].transform('first')
    td_data['flag_min_match_'+name2] = np.where(td_data['min_match_'+name2]==td_data['transcript_1'],1,0)
    td_data['flag_recip_min_match'] = np.where((td_data['flag_min_match_'+name1]==1)&
                                             (td_data['flag_min_match_'+name2]==1),1,0)
    td_data['num_recip_min_match_in_gene'] = td_data.groupby('gene_id')['flag_recip_min_match'].transform('sum')
    
    # Flag pairs that are ties in between at least 2 transcirpts in the other dataset
    td_data['flag_'+name1+'_tie'] = td_data[['transcript_1','prop_junction_diff','prop_ER_diff','prop_nt_diff','num_nt_diff']].duplicated(keep=False).astype(int)
    td_data['flag_'+name2+'_tie'] = td_data[['transcript_2','prop_junction_diff','prop_ER_diff','prop_nt_diff','num_nt_diff']].duplicated(keep=False).astype(int)

    # Flag FSM (all junctions matching) and ERM (all exonic regions shared) min matches
    td_data['flag_identical_recip_min_match'] = np.where((td_data['flag_recip_min_match']==1)&
                                                  (td_data['prop_nt_diff']==0),1,0)
    td_data['flag_FSM_recip_min_match'] = np.where((td_data['flag_recip_min_match']==1)&
                                            (td_data['prop_junction_diff']==0),1,0)
    td_data['flag_ERM_recip_min_match'] = np.where((td_data['flag_recip_min_match']==1)&
                                            (td_data['prop_ER_diff']==0),1,0)
    td_data['flag_ERM_noIR_recip_min_match'] = np.where((td_data['flag_recip_min_match']==1)&
                                            (td_data['prop_ER_diff']==0)&
                                            (td_data['num_IR_fragment_T1']+td_data['num_IR_fragment_T2']==0),1,0)
    td_data['flag_ERM_withIR_recip_min_match'] = np.where((td_data['flag_recip_min_match']==1)&
                                            (td_data['prop_ER_diff']==0)&
                                            (td_data['num_IR_fragment_T1']+td_data['num_IR_fragment_T2']>0),1,0)
    td_data['num_identical_recip_min_match_in_gene'] = td_data.groupby('gene_id')['flag_identical_recip_min_match'].transform('sum')
    td_data['num_FSM_recip_min_match_in_gene'] = td_data.groupby('gene_id')['flag_FSM_recip_min_match'].transform('sum')
    td_data['num_ERM_recip_min_match_in_gene'] = td_data.groupby('gene_id')['flag_ERM_recip_min_match'].transform('sum')
    td_data['num_ERM_noIR_recip_min_match_in_gene'] = td_data.groupby('gene_id')['flag_ERM_noIR_recip_min_match'].transform('sum')
    td_data['num_ERM_withIR_recip_min_match_in_gene'] = td_data.groupby('gene_id')['flag_ERM_withIR_recip_min_match'].transform('sum')
    td_data['prop_identical_recip_min_match_in_gene'] = td_data['num_identical_recip_min_match_in_gene']/td_data['num_recip_min_match_in_gene']
    td_data['prop_FSM_recip_min_match_in_gene'] = td_data['num_FSM_recip_min_match_in_gene']/td_data['num_recip_min_match_in_gene']
    td_data['prop_ERM_recip_min_match_in_gene'] = td_data['num_ERM_recip_min_match_in_gene']/td_data['num_recip_min_match_in_gene']
    td_data['prop_ERM_noIR_recip_min_match_in_gene'] = td_data['num_ERM_noIR_recip_min_match_in_gene']/td_data['num_recip_min_match_in_gene']
    td_data['prop_ERM_withIR_recip_min_match_in_gene'] = td_data['num_ERM_withIR_recip_min_match_in_gene']/td_data['num_recip_min_match_in_gene']
    
    # Flag different types of splicing differences in min matches
    #   (alt exon, alt donor/acceptors, IR)
    td_data['flag_nonoverlapping_recip_min_match'] = np.where((td_data['flag_recip_min_match']==1)&
                                                               (td_data['flag_nonoverlapping']==1))
    td_data['flag_alt_exon_recip_min_match'] = np.where((td_data['flag_recip_min_match']==1)&
                                                  (td_data['flag_alt_exon']==1),1,0)
    td_data['flag_alt_donor_acceptor_recip_min_match'] = np.where((td_data['flag_recip_min_match']==1)&
                                                  (td_data['flag_alt_donor_acceptor']==1),1,0)
    td_data['flag_IR_recip_min_match'] = np.where((td_data['flag_recip_min_match']==1)&
                                                  (td_data['flag_IR']==1),1,0)
    td_data['flag_5_variation_recip_min_match'] = np.where((td_data['flag_recip_min_match']==1)&
                                                             (td_data['flag_5_variation']==1),1,0)
    td_data['flag_3_variation_recip_min_match'] = np.where((td_data['flag_recip_min_match']==1)&
                                                             (td_data['flag_3_variation']==1),1,0)
    # Get reciprocal minimum match gene category
    conditionsRecip = [(td_data['transcript_in_gene']==name1+"_greater")&(td_data['num_recip_min_match_in_gene']==0),
                       (td_data['transcript_in_gene']==name1+"_greater")&(td_data['num_recip_min_match_in_gene']>0)&(td_data['num_recip_min_match_in_gene']<td_data['num_transcript_in_gene_'+name2]),
                       (td_data['transcript_in_gene']==name1+"_greater")&(td_data['num_recip_min_match_in_gene']>0)&(td_data['num_recip_min_match_in_gene']==td_data['num_transcript_in_gene_'+name2]),
                       (td_data['transcript_in_gene']==name2+"_greater")&(td_data['num_recip_min_match_in_gene']==0),
                       (td_data['transcript_in_gene']==name2+"_greater")&(td_data['num_recip_min_match_in_gene']>0)&(td_data['num_recip_min_match_in_gene']<td_data['num_transcript_in_gene_'+name1]),
                       (td_data['transcript_in_gene']==name2+"_greater")&(td_data['num_recip_min_match_in_gene']>0)&(td_data['num_recip_min_match_in_gene']==td_data['num_transcript_in_gene_'+name1]),
                       (td_data['transcript_in_gene']=="equal")&(td_data['num_recip_min_match_in_gene']==0),
                       (td_data['transcript_in_gene']=="equal")&(td_data['num_recip_min_match_in_gene']>0)&(td_data['num_recip_min_match_in_gene']<td_data['num_transcript_in_gene_'+name1]),
                       (td_data['transcript_in_gene']=="equal")&(td_data['num_recip_min_match_in_gene']>0)&(td_data['num_recip_min_match_in_gene']==td_data['num_transcript_in_gene_'+name1])]
    choicesRecrip = ["no_subset","partial_subset","subset","no_subset","partial_subset","subset","no_match","partial_match","match"]
    td_data['recip_min_pair_in_gene'] = np.select(conditionsRecip,choicesRecrip,"missing")
    if len(td_data[td_data['recip_min_pair_in_gene']=="missing"]) > 0:
        logger.error("Unexpected variable assignment for gene reciprocal minimum pair comparison")

    # Return minimum distance of transcript all pairs
    if all_pairs:
        return td_data
    # Return minimum distance of only minimum pairs of each transcript
    else:
        return td_data[(td_data['flag_min_match_'+name1]==1)|(td_data['flag_min_match_'+name2]==1)]