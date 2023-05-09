#!/usr/bin/env python

import argparse
import pandas as pd
import numpy as np
import os

import time

from trand.io import *

def getOptions():
        # Parse command line arguments
        parser = argparse.ArgumentParser(description="Consolidate GTF by identifying unique junction chains.")
        
        # Input data
        parser.add_argument("-g", "--gtf", dest="inGTF", required=True, help="Input GTF file.")
        parser.add_argument("-c", "--consol-prefix", dest="consolPrefix", required=False, default="tr", help="Prefix to use for consolidation transcript id values.")
        parser.add_argument("-v", "--verbose", dest="verbose", required=False, action="store_true", help="Verbose mode.")
        
        # Output data
        parser.add_argument("-o", "--output-prefix", dest="outPrefix", required=True, help="Full path and prefix for all output files.")
        
        args = parser.parse_args()
        return args


def create_junction_string(exon_data, verbose):
    """Create a junction strings for each transcript"""

    # Check for genes with inconsistent chromsome or strand
    exon_data = check_strand_and_chromosome(exon_data)
    
    transcript_groups = exon_data.groupby("transcript_id")

    # Iterate over transcript groups
    junction_chains = {}
    # num = 0
    for transcript_id, group in transcript_groups:
        # num = num + 1
        # if num > 100:
        #     break
        if len(group) > 1:
            # Sort exons by start position
            sorted_exons = group.sort_values(by="start").reset_index(drop=True)

            # Create junction strings using vectorized operations
            junctions = sorted_exons["seqname"].shift(-1).str.cat(
                    sorted_exons["end"].astype(int).astype(str), sep=":", na_rep=""
                ).str.cat(
                    sorted_exons["start"].shift(-1).fillna(0).astype(int).astype(str), sep=":", na_rep=""
                ).str.cat(
                    sorted_exons["strand"], sep=":", na_rep=""
                ).tolist()[:-1]

            # Concatenate junctions into a single string
            junction_chain = '|'.join(junctions)
            start = sorted_exons["start"].min()
            end = sorted_exons["end"].max()
        else:
            junction_chain = ""
            start = group["start"].iloc[0]
            end = group["end"].iloc[0]
        
        strand = group["strand"].iloc[0]
        seqname = group["seqname"].iloc[0]
        gene_id = group["gene_id"].iloc[0]
        junction_chains[transcript_id] = [junction_chain, transcript_id, gene_id, seqname, start, end, strand]
    return junction_chains

# all good, takes <10s on mel2mel (3million exons)
def check_strand_and_chromosome(exon_data):
        gene_groups = exon_data.groupby("gene_id")
        strand_check = gene_groups["strand"].nunique()
        
        if (strand_check > 1).any():
            bad_strand = list(strand_check[strand_check>1].index)
            for gene in bad_strand:
                print("!!! WARNING: gene {} contains transcripts/exons on both strands - removing from consolidation".format(gene))
            exon_data = exon_data[~exon_data["gene_id"].isin(bad_strand)]
        chr_check = gene_groups["seqname"].nunique()
        if (chr_check > 1).any():
            bad_chr = list(chr_check[chr_check>1].index)
            for gene in bad_chr:
                print("!!! WARNING: gene {} contains transcripts/exons on difference chromosomes - removing from consolidation".format(gene))
            exon_data = exon_data[~exon_data["gene_id"].isin(bad_chr)]
        return exon_data

def split_junction_string(df):
    # split junction string into a list of junctions
    junctions = df.set_index(["gene_id", "start", "end", "consolidation_transcript_id"])['junction_string'].str.split('|')
    
    # create empty arrays for exons and lengths
    exons = np.empty((0, 6), dtype=object)
    
    # iterate over junctions and split each into donor and acceptor coordinates
    for i in range(len(df)):
        # Check for monoexon transcripts (empty junction string)
        if type(junctions.iloc[i]) is not list:
            # Add exon to arrays
            exons = np.vstack([
                exons,
                np.array([df.loc[df.index[i], "seqname"],
                          df.loc[df.index[i], "start"],
                          df.loc[df.index[i], "end"],
                          df.loc[df.index[i], "strand"],
                          df.loc[df.index[i], "consolidation_transcript_id"],
                          df.loc[df.index[i], "gene_id"]])])
        else:
            # Get exon donor and acceptor values from junction coordinates
            junction_coords = np.array(junctions.iloc[i], dtype=object)
            donor_coords = np.array([val.split(':')[1] for val in junction_coords], dtype=int)
            acceptor_coords = np.array([val.split(':')[2] for val in junction_coords], dtype=int)
            
            # Calculate exon starts and ends
            exon_starts = np.concatenate(([df.loc[df.index[i], "start"]], acceptor_coords))
            exon_ends = np.concatenate((donor_coords, [df.loc[df.index[i], "end"]]))
            
            # Add exons to arrays
            exons = np.vstack([
                exons,
                np.concatenate(
                    (np.array([df.loc[df.index[i], "seqname"]]*len(exon_starts), dtype=object).reshape(-1,1),
                     exon_starts.reshape(-1,1),
                     exon_ends.reshape(-1,1),
                     np.array([df.loc[df.index[i], "strand"]]*len(exon_starts), dtype=object).reshape(-1,1),
                     np.array([df.loc[df.index[i], "consolidation_transcript_id"]]*len(exon_starts), dtype=object).reshape(-1,1),
                     np.array([df.loc[df.index[i], "gene_id"]]*len(exon_starts), dtype=object).reshape(-1,1),
                     ), axis=1)
                ])

    # Create dataframe of exons
    exons_df = pd.DataFrame(exons, columns=["seqname", "start", "end", "strand", "transcript_id", "gene_id"])
    return exons_df

def split_column_by_sep(df,col_name=None,sep=None,sort_list=None):
    # Split variable by some character like '|' or ',' and keep all other values the same
    if col_name == None:
        col_name = 'transcript_id'
    if sep == None:
        sep = "|"
    splitList = df[col_name].str.split(sep).apply(pd.Series, 1).stack()
    splitList.index = splitList.index.droplevel(-1)
    tempDF = df.copy()
    del(tempDF[col_name])
    splitDF = tempDF.join(splitList.rename(col_name))
    if sort_list != None:
        splitDF = splitDF.sort_values(by=sort_list)
    del(tempDF, splitList)
    return splitDF


def main():
    # Set consolidation prefix
    prefix = args.consolPrefix

    # Set up output file names
    outFiles = {
        "gtf": args.outPrefix + "_consolidated.gtf",
        "key": args.outPrefix + "_consolidation_key.csv"
    }

    # Get GTF junction strings
    ujcDict = create_junction_string(
        read_exon_data_from_file(args.inGTF),
        args.verbose
    )
    # ujcDf = create_junction_string(
    #     read_exon_data_from_file("/Volumes/blue/mcintyre/share/transcript_distance/wtc11_analysis/baseline_selection/PacBio_corrected_associated_gene_baseline_selection.gtf")
    # )
    # ujcDf = create_junction_string(
    #     read_exon_data_from_file("/Users/adalena/mclab/SHARE/McIntyre_Lab/useful_dmel_data/flybase617/flybase_files/dmel-all-r6.17.gtf")
    # )

    # Get monoexon transcript models
    # Create new dct from ujDct, copying all rows where the first value (junction string) is empty
    monoExon = {t:ujcDict[t] for t in ujcDict if ujcDict[t][0] == ""}
    if len(monoExon) > 0:
#        if args.verbose:
#            print("There are {} monoexon transcript")
        # Get non-overlapping genomic regions for monoexon transcripts
        monoexon_transcripts = pd.DataFrame(monoExon, index=pd.Index(["junction_string", "transcript_id", "gene_id", "seqname", "start", "end", "strand"])
                                            ).T.sort_values(by=["start", "end"])
        overlap = (monoexon_transcripts["start"].shift(-1) < monoexon_transcripts["end"]).cumsum()
        monoexon_transcripts["junction_string"] = overlap + 1

        # Get unique non-overlapping monoexon transcripts
        monoUJC = monoexon_transcripts.groupby(["gene_id","junction_string"]).agg({
            "seqname": "first",
            "start": "min",
            "end": "max",
            "strand": "first",
            "transcript_id": lambda x: "|".join(x)}).reset_index()
        del(monoexon_transcripts)

    # Get multiexon transcript models
    multiExon = {t:ujcDict[t] for t in ujcDict if ujcDict[t][0] != ""}
    del(ujcDict)

    if len(multiExon) > 0:
        # Get unique junction chains of multiexon transcript models
        multiUJC = pd.DataFrame(multiExon, index=pd.Index(["junction_string", "transcript_id", "gene_id", "seqname", "start", "end", "strand"])
                            ).T.groupby(["gene_id","junction_string"]).agg({
                                "seqname": "first",
                                "start": "min",
                                "end": "max",
                                "strand": "first",
                                "transcript_id": lambda x: "|".join(x)}).reset_index()

    # Merge monoexon and multiexon transcript model information
    if len(monoExon) > 0 and len(multiExon) > 0:
        allUJC = pd.concat([monoUJC, multiUJC], ignore_index=True)
        del(monoUJC)
        del(multiUJC)
    elif len(monoExon) > 0 and len(multiExon) == 0:
        allUJC = monoUJC.copy()
        del(monoUJC)
    else:
        allUJC = multiUJC.copy()
        del(multiUJC)
                                
    allUJC["consol_transcript_length"] = allUJC["end"] - allUJC["start"]
#    allUJC.sort_values(by=["gene_id", "consol_transcript_length", "transcript_id"], ascending=[True, False, True], inplace=True)
#    allUJC.sort_values(by=["gene_id", "consol_transcript_length"], ascending=[True, False], inplace=True)
#    allUJC.sort_values(by=['gene_id', 'consol_transcript_length', 'start', 'junction_string'], ascending=[False, False, True, True], inplace=True)

    sort_order = {"gene_id": "asc", "consol_transcript_length": "asc", "start": "asc", "transcript_id": "asc"}
    allUJC = allUJC.sort_values(by=list(sort_order.keys()), ascending=[True if val=="asc" else False for val in sort_order.values()])

    allUJC["transcript_rank_in_gene"] = (
        allUJC.groupby("gene_id")["consol_transcript_length"].rank(method="first")
    )
    allUJC["consolidation_transcript_id"] = (
        prefix
        + "_"
        + allUJC["gene_id"]
        + "_"
        + allUJC["transcript_rank_in_gene"].astype(int).map(str)
    )

    # Make output GTF file of consolidated models
    if os.path.isfile(outFiles["gtf"]):
        os.remove(outFiles["gtf"])
    write_gtf(
        split_junction_string(allUJC),
        outFiles,
        "gtf"
    )

    # Make key file
    keys = split_column_by_sep(allUJC[["gene_id", "transcript_id", "consolidation_transcript_id"]], col_name="transcript_id", sep="|")
    keys["num_transcript_in_consol_transcript"] = keys.groupby(
            "consolidation_transcript_id")["transcript_id"].transform('nunique')
    keys["num_transcript_id_in_gene"] = keys.groupby(
            "gene_id")["transcript_id"].transform('nunique')
    keys["num_consol_transcript_id_in_gene"] = keys.groupby(
            "gene_id")["consolidation_transcript_id"].transform('nunique')
    keys["flag_gene_consolidated"] = np.where(
        keys["num_transcript_id_in_gene"] > keys["num_consol_transcript_id_in_gene"],
        1,
        0
    )
    keys.to_csv(outFiles["key"], index=False)


if __name__ == '__main__':
    # Parse command line arguments
    global args
    
    tic = time.perf_counter()
    args = getOptions()
    main()
    toc = time.perf_counter()
    print(f"Complete! Operation took {toc-tic:0.4f} seconds.")
