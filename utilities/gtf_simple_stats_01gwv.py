# -*- coding: utf-8 -*-
"""
Created on Mon Feb 13 15:10:58 2022

@author: gvanveckhoven
"""

import argparse
import pandas as pd
import os

#Import trand functions
from trand import io
from trand.event_analysis import nCr

def getOptions():
    # Parse command line arguments
    parser = argparse.ArgumentParser(
        description=(
                "Input either one GTF file or two GTF files to find the pairs between the genes of both or wihtin the genes of one"
        )
    )
    
    #input data
    parser.add_argument( 
        "-g1",
        "--gtf1",
        dest="gtf1",
        required=True,
        help="Input one GTF file to calculate the pairs of transcripts in the one," 
             "otherwise it will be used to find the pairs between the first inputted" 
             "GTF file and the second"
    )
    parser.add_argument( 
        "-g2",
        "--gtf2",
        dest="gtf2",
        help="Second GTF file to calculate pairs with the first"
    )
    
    # Output arguments
    parser.add_argument(
        "-o", 
        "--output", 
        dest="outFile", 
        default= (os.getcwd()+'/stats_info.csv'),
        help="Output file name for subset GTF"
    )
    
    args = parser.parse_args()
    return args

def validate_input(gtf):
    #takes GTF file inputs and creates df's
    gtf_df = io.read_exon_data_from_file(gtf)
    return gtf_df
    


def pairs1(gtf1_df):
    #finds the number of pairs in the GTF dataframe
    gtf1_df = gtf1_df.groupby("gene_id")["transcript_id"].nunique().reset_index()
    gtf1_df = gtf1_df.drop(gtf1_df[gtf1_df['transcript_id'] < 2].index)
    total_pairs = gtf1_df['transcript_id'].map(lambda x: nCr (x,2)).sum()
    return (int(total_pairs))
 
def pairs2(gtf1_df, gtf2_df):
    #finds the number of pairs between two GTF dataframes
    gtf1_df = gtf1_df.groupby("gene_id")["transcript_id"].nunique().reset_index()
    gtf2_df = gtf2_df.groupby("gene_id")["transcript_id"].nunique().reset_index()
    
    gtf1_df['roz'] = gtf1_df['gene_id'].map(gtf2_df.set_index('gene_id')['transcript_id'])
    gtf1_df['pairs'] = gtf1_df['roz'] * gtf1_df['transcript_id']
    total_pairs_between = gtf1_df['pairs'].sum()
    return (total_pairs_between)

def main():
    #Convert gtf to df and stats calculated
    gtf1_df = validate_input(args.gtf1)
    gtf1_genes = gtf1_df.loc[:, 'gene_id'].nunique()
    gtf1_transcripts = gtf1_df.loc[:, 'transcript_id'].nunique()
    gtf1_chromosomes = gtf1_df.loc[:, 'seqname'].nunique()
    gtf1_pairs = pairs1(gtf1_df)

    if args.gtf1 is not None and args.gtf2 is None:
        #Iff one gtf file is input, the stats are added to output csv
        single_info = [gtf1_genes, gtf1_transcripts, gtf1_chromosomes, gtf1_pairs]
        single_info_df = pd.DataFrame(data = single_info, index = ['genes', 'transcripts', 'chromosomes', 'pairs'])
        single_info_df.to_csv(args.outFile, header=False) 

    if args.gtf2 is not None:
        #if a second gtf file is inputted then a different set of values is added to the output csv
        gtf2_df = validate_input(args.gtf2)
        gtf2_genes = gtf2_df.loc[:, 'gene_id'].nunique()
        gtf2_transcripts = gtf2_df.loc[:, 'transcript_id'].nunique()
        gtf2_chromosomes = gtf2_df.loc[:, 'seqname'].nunique()
        gtf2_pairs = pairs1(gtf2_df)
        
    if args.gtf1 and args.gtf2 is not None:
        gtf_pairs_between = pairs2(gtf1_df, gtf2_df)
        dual_info = [gtf1_genes, gtf2_genes, gtf1_transcripts, gtf2_transcripts, 
                     gtf1_chromosomes, gtf2_chromosomes, gtf1_pairs, gtf2_pairs, 
                     gtf_pairs_between]
        dual_info_df = pd.DataFrame(data = dual_info, index = ['genes1', 'genes2', 'transcripts1', 'transcripts2',
                                                'chromsomes1', 'chromosomes2', 'pairs1', 'pairs2',
                                                'pairs_between'])
        dual_info_df.to_csv(args.outFile, header=False)    

if __name__ == '__main__':
    # Parse command line arguments
    global args
    args = getOptions()
    main()
