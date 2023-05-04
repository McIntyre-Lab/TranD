# -*- coding: utf-8 -*-
"""
Created on Mon Feb 13 15:10:58 2022

@author: gvanveckhoven
"""

import argparse
import pandas as pd
import os
from collections import Counter

#Import trand functions
import trand
from trand import io
from trand import event_analysis
from trand.event_analysis import nCr

def getOptions():
        """
        
        Function to store user input via argparse
        
        Returns
        -------
        args : ARGPARSE ARGUMENTS
                User input via argparse.
        
        """
        
        # Parse command line arguments
        parser = argparse.ArgumentParser(
            description=(
                    "Input either 1 GTF file (--gtf1) or add another optional 2nd GTF file (--gtf2) with the same geneIDs."
                    "Will output the following file: a table summarizing the number of genes, number of transcripts, "
                    "number of chromosomes and number of transcript pairs (within 1 GTF for across 2 GTF files)."
                    "Also has the option to output 2 lists of genes based on the number of transcripts per gene "
                    "specified by the user (--cutoff). One file >= to cutoff and second file < cutoff."
            )
        )
        
        #input data
        parser.add_argument( 
                "-g1",
                "--gtf1",
                dest="gtf1",
                required=True,
                help="Input one GTF file to calculate the pairs of transcripts by gene."
        )
        
        parser.add_argument( 
                "-g2",
                "--gtf2",
                dest="gtf2",
                help="Input a second GTF file to calculate pairs of transcripts across the 2 GTF files."
        )
        
        parser.add_argument(
                "-cut",
                "--cutoff",
                dest="cutoff",
                help="Input a cutoff value for the number of transcripts per gene. Used to generate 2 lists of genes for each GTF "
                      "based on the number of transcripts per gene." 
        )
        
        # Output arguments
        parser.add_argument(
            "-o", 
            "--output", 
            dest="outFile", 
            default= (os.getcwd()+'/stats_info.csv'),
            help="Enter an filename with apath for the output summary table."
            "Defaults to the parent directory with the name stats_info.csv."
        )
        
        args = parser.parse_args()
        return args

def validateInput(gtf):
        """
        
        Takes GTF file inputs and creates dataframes
        
        Parameters
        ----------
        gtf : FILE
                GTF file input.
        
        Returns
        -------
        gtf_df : DATAFRAME
                A dataframe that organizes the information in the GTF in a usable way.
        """
        #takes GTF file inputs and creates df's
        gtf_df = io.read_exon_data_from_file(gtf)
        return gtf_df
    


def calc1GTFPair(gtf_df):
        """
        
        Finds the number of transcript pairs in the dataframe.
        
        Parameters
        ----------
        gtf1_df : DATAFRAME
                A dataframe created from a GTF file.
        
        Returns
        -------
        INT
                Total number of pairs of transcripts in the GTF file.
        
        """
        #finds the number of pairs in the GTF dataframe
        gtf_df = gtf_df.groupby("gene_id")["transcript_id"].nunique().reset_index()
        gtf_df = gtf_df.drop(gtf_df[gtf_df['transcript_id'] < 2].index)
        total_pairs = gtf_df['transcript_id'].map(lambda x: nCr (x,2)).sum()
        return (int(total_pairs))
         
def calc2GTFPair(gtf1_df, gtf2_df):
        """
        
        Finds the number of transcript pairs between both dataframes.

        Parameters
        ----------
        gtf1_df : DATAFRAME
                A dataframe created from a GTF file.
        gtf2_df : DATAFRAME
                Another dataframe created from a GTF file.

        Returns
        -------
        total_pairs_between : INT
                Total number of pairs of transcripts between both GTF files.

        """
        #finds the number of pairs between two GTF dataframes
        gtf1_df = gtf1_df.groupby("gene_id")["transcript_id"].nunique().reset_index()
        gtf2_df = gtf2_df.groupby("gene_id")["transcript_id"].nunique().reset_index()
        
        gtf1_df['roz'] = gtf1_df['gene_id'].map(gtf2_df.set_index('gene_id')['transcript_id'])
        gtf1_df['pairs'] = gtf1_df['roz'] * gtf1_df['transcript_id']
        total_pairs_between = gtf1_df['pairs'].sum()
        return (total_pairs_between)

def sortGenesByCutoff(gtf1_df, gtf2_df, cutoff):
        """
        
        Outputs two lists of genes that have above and below a certain number of
        transcripts.
        
        Works for one or two GTFs.

        Parameters
        ----------
        gtf1_df : DATAFRAME
                A dataframe created from a GTF file.
        gtf2_df : DATAFRAME
                Another dataframe created from a GTF file.
        cutoff : INT
                The transcript threshold.

        Returns
        -------
        aboveDf : DATAFRAME
                A list of genes with a number of transcripts above the cutoff.
        belowDf : DATAFRAME
                A list of genes with a number of transcripts below the cutoff.

        """
        if gtf2_df is not None:
                tempDf1 = gtf1_df.drop_duplicates(subset='transcript_id')
                tempDf2 = gtf2_df.drop_duplicates(subset='transcript_id')
                
                geneDct = Counter(tempDf1['gene_id']) + Counter(tempDf2['gene_id'])

        else:
                tempDf = gtf1_df.drop_duplicates(subset='transcript_id')
                geneDct = Counter(tempDf['gene_id'])
        
        aboveLst = []
        belowLst = []
        
        for key, value in geneDct.items():
                
                if (value >= cutoff):
                        aboveLst.append(key)
                else:
                        belowLst.append(key)
        
        aboveDf = pd.DataFrame(aboveLst)
        belowDf = pd.DataFrame(belowLst)

        return aboveDf, belowDf

def main():
        """
        
        Run the program. (calls the above functions)

        Returns
        -------
        Nothing.

        """
        #Convert gtf to df and stats calculated
        gtf1_df = validateInput(gtf=args.gtf1)
        gtf1_genes = gtf1_df.loc[:, 'gene_id'].nunique()
        gtf1_transcripts = gtf1_df.loc[:, 'transcript_id'].nunique()
        gtf1_chromosomes = gtf1_df.loc[:, 'seqname'].nunique()
        gtf1_pairs = calc1GTFPair(gtf_df=gtf1_df)
        
        if args.gtf1 is not None and args.gtf2 is None:
            #Iff one gtf file is input, the stats are added to output csv
            single_info = [gtf1_genes, gtf1_transcripts, gtf1_chromosomes, gtf1_pairs]
            single_info_df = pd.DataFrame(data = single_info, index = ['genes', 'transcripts', 'chromosomes', 'pairs'])
            single_info_df.to_csv(args.outFile, header=False) 
        
        gtf2_df = None
        if args.gtf2 is not None:
            #if a second gtf file is inputted then a different set of values is added to the output csv
            gtf2_df = validateInput(gtf=args.gtf2)
            gtf2_genes = gtf2_df.loc[:, 'gene_id'].nunique()
            gtf2_transcripts = gtf2_df.loc[:, 'transcript_id'].nunique()
            gtf2_chromosomes = gtf2_df.loc[:, 'seqname'].nunique()
            gtf2_pairs = calc1GTFPair(gtf_df=gtf2_df)
            
        if args.gtf1 and args.gtf2 is not None:
            gtf_pairs_between = calc2GTFPair(gtf1_df=gtf1_df, gtf2_df=gtf2_df)
            dual_info = [gtf1_genes, gtf2_genes, gtf1_transcripts, gtf2_transcripts, 
                         gtf1_chromosomes, gtf2_chromosomes, gtf1_pairs, gtf2_pairs, 
                         gtf_pairs_between]
            
            dual_info_df = pd.DataFrame(data = dual_info, index = ['genes1', 'genes2', 'transcripts1', 'transcripts2',
                                                    'chromsomes1', 'chromosomes2', 'pairs1', 'pairs2',
                                                    'pairs_between'])
            
            dual_info_df.to_csv(args.outFile, header=False)
        
        if args.cutoff is not None:
                aboveDf, belowDf = sortGenesByCutoff(gtf1_df=gtf1_df, gtf2_df=gtf2_df, cutoff=int(args.cutoff))
                
                filepath = os.path.dirname(args.outFile)
                
                aboveDf.to_csv(filepath + "/genes_above_" + str(args.cutoff) + ".csv", index=False, header=False)
                belowDf.to_csv(filepath + "/genes_below_" + str(args.cutoff) + ".csv", index=False, header=False)


if __name__ == '__main__':
    # Parse command line arguments
    global args
    args = getOptions()
    main()
    print("Complete!")
