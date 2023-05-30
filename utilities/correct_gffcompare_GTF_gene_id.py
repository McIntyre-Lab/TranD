#!/usr/bin/env python

import argparse
import pandas as pd
import numpy as np
import csv
import sys


def getOptions():
    # Parse command line arguments
    parser = argparse.ArgumentParser(description="Correct the gene_id values of the GTF-like output from GFFcompare (v12.2 or higher) to match the gene_id values associated with the transcripts")

    # Input data
    parser.add_argument(
        "-a",
        "--annotation",
        dest="inAnnot",
        required=True,
        help=(
            "Input reference annotated GTF-like file from GFFcompare >= v12.2 (*.annotated.gtf). "
            "NOTE: If multiple GTF files need to have gene_id values assigned, must provide "
            "gffcompare with a single concatenated GTF where each individual file has a sampleID "
            "in the \'source\' column (column 2) of each GTF file."
        )
    )
    parser.add_argument(
        "-k",
        "--keep-names",
        dest="inKeep",
        required=False,
        action="store_true",
        help=(
            "When a reference gene_id was not associated, use original gene_id "
            "(Default: Replace with GFFcompare XLOC values)."
        )
    )
    parser.add_argument(
        "-n",
        "--use-gene_name",
        dest="inName",
        required=False,
        action="store_true",
        help=(
            "When a ref_gene_id is not present in GFFcompare output for class "
            "codes not in [s, x, y, i, p, r, u], check for ref_gene_name or "
            "gene_name to use instead (Default: Only use ref_gene_id, else "
            "replace with GFFcompare XLOC values or original gene_id if -k is used)."
        )
    )
    parser.add_argument(
        "-g",
        "--gtf",
        dest="inGTF",
        required=True,
        help=(
            "Input GTF file to be modified (original GTF given to GFFcompare). "
            "NOTE: The name of this file MUST match what was given to "
            "gffcompare so the correct transcripts can be identified in the "
            "output if multiple files were provided to gffcompare."
        )
    )

    # Output data
    parser.add_argument(
        "-o",
        "--output",
        dest="outFile",
        required=True,
        help="Output GTF file name for isoform exons with associated gene_id (file is UNSORTED)"
    )
    parser.add_argument(
        "--key",
        dest="outKey",
        required=True,
        help=(
            "Output CSV key file for each transcript_id, input_gene_id, and "
            "output_gene_id, where transcript_id is the list of all unique "
            "transcripts in the output GTF file, input_gene_id is the gene_id "
            "of the transcript in the input GTF file, and output_gene_id is the "
            "new associated gene_id based on the gffcompare annotation file and "
            "associated with the transcript_id in the output GTF file."
        )
    )

    args = parser.parse_args()
    return args

def main():

    # Get gffcompare annotation output 
    annotDF = pd.read_csv(args.inAnnot,names=['chr','source','feature','start','end','score',
                                        'strand','frame','attribute'], sep="\t",low_memory=False)
    print("Total lines in gffcompare output = "+str(len(annotDF)))
    print("Total transcripts in gffcompare output = "+str(len(annotDF[annotDF["feature"]=="transcript"])))
    print("Total input sources contributing to gffcompare output = "+str(annotDF["source"].nunique()))

    # Get gene_id and transcript_id values from attributes column
    # Get transcript_id to gene_id pairs using transcript features of gffcompare output
    # If ref_gene_id attribute present (from gffcompare v12.2), assign as gene_id
    # else if xloc attribute present, assign as gene_id (if requested)
    # else not properly formatted
    xcrptDF = annotDF[annotDF["feature"]=="transcript"].reset_index()
    for i in xcrptDF.index:
        raw_attrs = xcrptDF.at[i, 'attribute']
        attr_list = [x.strip() for x in raw_attrs.strip().split(';')]
        ## previous gffcompare versions do not have "ref_gene"
        exp_attrs = [x for x in attr_list if 'transcript_id' in x or 'gene_id' in x or 'xloc' in x or "gene_name" in x or "class_code" in x]
        transcript_id, orig_gene_id, ref_gene_id, ref_gene_name, gene_name, xloc, classCode = np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan
        for item in exp_attrs:
            if 'transcript_id' in item:
                transcript_id = item.split('transcript_id')[-1].strip().strip('\"')
            elif 'ref_gene_id' in item:
                ref_gene_id = item.split('ref_gene_id')[-1].strip().strip('\"')
            elif 'ref_gene_name' in item:
                ref_gene_name = item.split('ref_gene_name')[-1].strip().strip('\"')
            elif 'gene_name' in item:
                gene_name = item.split('gene_name')[-1].strip().strip('\"')
            elif 'gene_id' in item:
                orig_gene_id = item.split('gene_id')[-1].strip().strip('\"')
            elif 'xloc' in item:
                xloc = item.split('xloc')[-1].strip().strip('\"')
            elif 'class_code' in item:
                classCode = item.split('class_code')[-1].strip().strip('\"')
        if transcript_id is np.nan:
            print("WARNING: transcript_id not found in {}".format(xcrptDF[i]))
        keepOrigGene = args.inKeep
        useGeneName = args.inName
        if ref_gene_id is np.nan:
            if not useGeneName:
                if not keepOrigGene:
                    if xloc is np.nan:
                        print("WARNING: no attribute for gene_id assignment (ref_gene_id or xloc) found in {}".format(xcrptDF[i]))
                    else:
                        gene_id = xloc
                else:
                    gene_id = orig_gene_id
            else:
                if classCode not in ['s', 'x', 'y', 'i', 'p', 'r', 'u']:
                    if ref_gene_name is np.nan and gene_name is np.nan:
                        if not keepOrigGene:
                            if xloc is np.nan:
                                print("WARNING: no attribute for gene_id assignment (ref_gene_id or xloc) found in {}".format(xcrptDF[i]))
                            else:
                                gene_id = xloc
                        else:
                            gene_id = orig_gene_id
                    else:
                        if ref_gene_name is np.nan:
                            gene_id = gene_name
                        else:
                            gene_id = ref_gene_name
                else:
                    if not keepOrigGene:
                        if xloc is np.nan:
                            print("WARNING: no attribute for gene_id assignment (ref_gene_id or xloc) found in {}".format(xcrptDF[i]))
                        else:
                            gene_id = xloc
                    else:
                        gene_id = orig_gene_id
        else:
            gene_id = ref_gene_id
        xcrptDF.at[i, 'gene_id'] = str(gene_id)
        xcrptDF.at[i, 'transcript_id'] = str(transcript_id)

    # Get GTF file and count total lines
    gtf = pd.read_csv(args.inGTF,names=['chr','source','feature','start','end','score',
                                        'strand','frame','attribute'], sep="\t",low_memory=False,comment="#")
    print("Total lines in original GTF= "+str(len(gtf)))

    # Check for unique pairs of gene and transcript
    # Duplicate values could be due to trans-spliced inputs
    if xcrptDF["transcript_id"].duplicated().any():
        print("!!!ERROR: Duplicate gene_id:transcript_id pairs in annotation file"
              "...can be due to trans-spliced inputs.\n"
              "\tThe following are the rows with duplicated gene-transcript:\n{}".format(
                  xcrptDF[xcrptDF["transcript_id"].duplicated(keep=False)].to_string()))
        exit()

    # Get gene_id and transcript_id values from attributes column
    for i in gtf.index:
        raw_attrs = gtf.at[i, 'attribute']
        attr_list = [x.strip() for x in raw_attrs.strip().split(';')]
        g_t_attrs = [x for x in attr_list if 'transcript_id' in x or 'gene_id' in x]
        gene_id, transcript_id = np.nan, np.nan
        for item in g_t_attrs:
            if 'gene_id' in item:
                gene_id = item.split('gene_id')[1].strip().strip('\"')
                gene_idx = g_t_attrs.index(item)
            elif 'transcript_id' in item:
                transcript_id = item.split('transcript_id')[-1].strip().strip('\"')
                transcript_idx = g_t_attrs.index(item)
        if transcript_id == np.nan and gtf.at[i, 'feature'] != "gene" :
            print("WARNING: transcript_id not found in {}".format(gtf[i]))
        gtf.at[i, 'gene_id'] = str(gene_id)
        gtf.at[i, 'gene_idx'] =  int(gene_idx)
        gtf.at[i, 'transcript_id'] = str(transcript_id)
        gtf.at[i, 'transcript_idx'] = int(transcript_idx)
    print("Total transcripts in original GTF= "+str(gtf["transcript_id"].nunique()))
    print("Total genes in original GTF= "+str(gtf["gene_id"].nunique()))
    
    # Merge gene_id from gffcompare output to original gtf on transcript_id
    mergeDFxcrpt = pd.merge(
            gtf,
            xcrptDF[["transcript_id", "gene_id"]],
            how="left",
            on="transcript_id",
            suffixes=["_gtf", "_cmp"],
            indicator="merge_check",
            validate="m:1"
    )
    if mergeDFxcrpt["merge_check"].value_counts()["left_only"] > 0:
        print("WARNING: Original GTF contains transcripts not present in gffcompare output.")
    
    # Check for proper merge
    if len(mergeDFxcrpt[mergeDFxcrpt['gene_id_cmp'].isna()]) > 0:
        sys.exit("!!! ERROR : UNEXPECTED MERGE RESULTING IN MISSING gene_id ASSOCIATION")
    
    # Create new attribute column with associated_gene as gene_id
    for i in mergeDFxcrpt.index:
        if mergeDFxcrpt.at[i, "gene_idx"] == 0:
            new_attribute = "gene_id \"" + mergeDFxcrpt.at[i, "gene_id_cmp"] + "\";" + ";".join(mergeDFxcrpt.at[i, "attribute"].split(";")[1:])
        else:
            new_attribute = ";".join(mergeDFxcrpt.at[i, "attribute"].split(";")[0:int(mergeDFxcrpt.at[i, "gene_idx"])]) + "; gene_id \"" + mergeDFxcrpt.at[i, "gene_id_cmp"] + "\";" + ";".join(mergeDFxcrpt.at[i, "attribute"].split(";")[int(mergeDFxcrpt.at[i, "gene_idx"])+1:])
        mergeDFxcrpt.at[i, "new_attribute"] = new_attribute

    # Output GTF file with associated_gene in new attribute
    mergeDFxcrpt[['chr','source','feature','start','end','score','strand','frame',
             'new_attribute']].to_csv(args.outFile,sep="\t",index=False,header=False,
             doublequote=False,quoting=csv.QUOTE_NONE)

    # Output key file for transcript_id to input_gene_id to output_gene_id
    keyDF = mergeDFxcrpt[["transcript_id", "gene_id_gtf", "gene_id_cmp"]].drop_duplicates().rename(
        columns={"gene_id_gtf": "input_gene_id","gene_id_cmp": "output_gene_id"})
    keyDF.to_csv(args.outKey, index=False)


if __name__ == '__main__':
    # Parse command line arguments
    global args
    args = getOptions()
    main()

