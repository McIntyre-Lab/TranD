#!/usr/bin/env python

import argparse
import pandas as pd
import numpy as np
import os

def list_files(directory, extension, allow_ext=False):
    if allow_ext == True:
        return (f for f in os.listdir(directory) if f.endswith('.' + extension) or "."+extension+"." in f)
    else:
        return (f for f in os.listdir(directory) if f.endswith('.' + extension))

def getOptions():
    # Parse command line arguments
    parser = argparse.ArgumentParser(description="Describe complexity of transcriptome(s) for a give GTF or set of GTF files")

    # Input data
    parser.add_argument("-i", "--input-directory", dest="inDir", required=True, help="Input directory of GTF files (all *.gtf files in directory will be processed)")
    parser.add_argument("-f", "--force", dest="force", required=False, action='store_true', help="Force output to be created for all *.gtf files in the directory, without this any with output already detected will be skipped")
    parser.add_argument("-e", "--exon-per-gene", dest="inExonPerGene", required=False, action='store_true', help="Output the number of unique exons within each gene")

    # Output data
    parser.add_argument("-o", "--output-directory", dest="outDir", required=True, help="Output directory for counts")

    args = parser.parse_args()
    return args

def main():
    # Loop over *.gtf files in input directory
    for file in list_files(args.inDir, "gtf", allow_ext=True):
        # Check if output was already generated
        if os.path.basename(file).split(".")[-1] == "gtf":
            prefix = ".".join(os.path.basename(file).split(".")[:-1])
        elif os.path.basename(file).split(".")[-2] == "gtf":
            prefix = ".".join(os.path.basename(file).split(".")[:-2])
        else:
            print("For file {}...\n\tERROR: Unrecognized file extension, skipping file".format(file))
            continue
        if args.force is False and os.path.isfile("{}/{}.transcriptome_counts.csv".format(args.outDir,prefix)):
            print("For file {}...\n\tWARNING: Transcriptome counts output file detected, skipping GTF".format(file))
            continue
        
        # Get GTF file
        gtf = pd.read_csv(args.inDir+"/"+ file,names=['chr','source','feature',
                                                      'start','end','score',
                                                      'strand','frame','attribute'],dtype=str,compression='infer',sep="\t",low_memory=False)
        # Remove header/comment lines if present
        gtf = gtf[~gtf['chr'].str.startswith("#")]
        
        # Select only exon features
        gtf = gtf[gtf['feature']=='exon']

        # Get attribute labels 1 and 2 (name 3 is used if the gene name has a space in it)
        gtf['attribute_name_1'] = gtf['attribute'].str.split(" ").str[0]
        gtf['attribute_name_2'] = gtf['attribute'].str.split(" ").str[2]
        gtf['attribute_name_3'] = gtf['attribute'].str.split(" ").str[3]
        
        
        # Check that attributes are only gene_id and transcript_id
        if len(gtf['attribute_name_1'].unique()) > 2 or len(gtf['attribute_name_2'].unique()) > 2 :
            if len(gtf[(gtf['attribute_name_1']!="transcript_id")&(gtf['attribute_name_2']!="transcript_id")&(gtf['attribute_name_3']!="transcript_id")]) > 0:
                print("For file {}...\n\tERROR: First two attributes contain more than 2 types, skipping GTF".format(file))
        if 'transcript_id' not in gtf['attribute_name_1'].unique() and 'transcript_id' not in gtf['attribute_name_2'].unique():
            print("For file {}...\n\tERROR: transcript_id not contained within first two attributes, skipping GTF".format(file))
        if 'gene_id' not in gtf['attribute_name_1'].unique() and 'gene_id' not in gtf['attribute_name_2'].unique():
            print("For file {}...\n\tERROR: gene_id not contained within first two attributes, skipping GTF".format(file))

        # Extract transcript_id and gene_id from attribute column
        gtf['transcript_id'] = np.where(gtf['attribute_name_1']=="transcript_id",
                                           gtf['attribute'].str.split(";").str[0].str.split("\"").str[1],
                                           np.where((gtf['attribute_name_2']=="transcript_id")|(gtf['attribute_name_3']=="transcript_id"),
                                                    gtf['attribute'].str.split(";").str[1].str.split("\"").str[1],"oops"))
        gtf['gene_id'] = np.where(gtf['attribute_name_1']=="gene_id",
                                           gtf['attribute'].str.split(";").str[0].str.split("\"").str[1],
                                           np.where(gtf['attribute_name_2']=="gene_id",
                                                    gtf['attribute'].str.split(";").str[1].str.split("\"").str[1],"oops"))

        # Get transcript and gene level counts
        gtf['num_exon'] = 1
        gtf['exon_id'] = gtf['chr'].map(str)+":"+gtf['start'].map(str)+":"+gtf['end'].map(str)+":"+gtf['strand'].map(str)
        transcriptDF = gtf.groupby('transcript_id').agg({'num_exon':'sum'}).reset_index()
        geneDF = gtf.groupby('gene_id').agg({'transcript_id':'nunique','num_exon':'sum','exon_id':'nunique'}).reset_index().rename(columns={'transcript_id':'num_transcript','exon_id':'num_uniq_exon'})

        # Output gene level files if requested
        if args.inExonPerGene == True:
            geneDF[['gene_id','num_uniq_exon']].to_csv("{}/{}.uniqExonPerGene.csv".format(args.outDir,prefix),index=False)

        # Get number of transcripts
        counts = pd.DataFrame()
        counts.loc[0,'num_transcript'] = geneDF['num_transcript'].sum()
        
        # Get number of genes
        counts.loc[0,'num_gene'] = len(geneDF)
        
        # Get number of exons (same exon in multiple transcripts counted multiple times)
        counts.loc[0,'num_exon'] = geneDF['num_exon'].sum()
        
        # Get number of unique exons (can still overlap but not have the same coordinates)
        # Must use GTF and not gene level counts since unique exons can be in multiple genes
        counts.loc[0,'num_uniqExon'] = len(gtf[['chr','start','end','strand']].drop_duplicates())
        
        # Get minimum values
        counts.loc[0,'min_transcriptPerGene'] = geneDF['num_transcript'].min()
        counts.loc[0,'min_exonPerTranscript'] = transcriptDF['num_exon'].min()
        counts.loc[0,'min_exonPerGene'] = geneDF['num_uniq_exon'].min()
        
        # Get Q1 values
        counts.loc[0,'q1_transcriptPerGene'] = geneDF['num_transcript'].quantile(0.25)
        counts.loc[0,'q1_exonPerTranscript'] = transcriptDF['num_exon'].quantile(0.25)
        counts.loc[0,'q1_exonPerGene'] = geneDF['num_uniq_exon'].quantile(0.25)
        
        # Get median values
        counts.loc[0,'median_transcriptPerGene'] = geneDF['num_transcript'].median()
        counts.loc[0,'median_exonPerTranscript'] = transcriptDF['num_exon'].median()
        counts.loc[0,'median_exonPerGene'] = geneDF['num_uniq_exon'].median()

        # Get Q3 values
        counts.loc[0,'q3_transcriptPerGene'] = geneDF['num_transcript'].quantile(0.75)
        counts.loc[0,'q3_exonPerTranscript'] = transcriptDF['num_exon'].quantile(0.75)
        counts.loc[0,'q3_exonPerGene'] = geneDF['num_uniq_exon'].quantile(0.75)
        
        # Get max values
        counts.loc[0,'max_transcriptPerGene'] = geneDF['num_transcript'].max()
        counts.loc[0,'max_exonPerTranscript'] = transcriptDF['num_exon'].max()
        counts.loc[0,'max_exonPerGene'] = geneDF['num_uniq_exon'].max()
        
        # Get mean values
        counts.loc[0,'mean_transcriptPerGene'] = geneDF['num_transcript'].mean()
        counts.loc[0,'mean_exonPerTranscript'] = transcriptDF['num_exon'].mean()
        counts.loc[0,'mean_exonPerGene'] = geneDF['num_uniq_exon'].mean()

        # Get std
        counts.loc[0,'std_transcriptPerGene'] = geneDF['num_transcript'].std()
        counts.loc[0,'std_exonPerTranscript'] = transcriptDF['num_exon'].std()
        counts.loc[0,'std_exonPerGene'] = geneDF['num_uniq_exon'].std()
        
        # Output counts to file
        counts.to_csv("{}/{}.transcriptome_counts.csv".format(args.outDir,prefix),index=False)

if __name__ == '__main__':
    # Parse command line arguments
    global args
    args = getOptions()
    main()
