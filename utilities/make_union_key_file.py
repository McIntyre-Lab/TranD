#!/usr/bin/env python

import argparse
import pandas as pd

def getOptions():
    # Parse command line arguments
    parser = argparse.ArgumentParser(description="Make union key file using UJC of 1) dataset1, 2) "
                                     "dataset2, and 3) the concatenation of consolidated dataset1 "
                                     "and consolidated dataset2.")

    # Input data
    parser.add_argument(
        "-1",
        "--dataset1",
        dest="inFile1",
        required=True,
        help="UJC key file for dataset 1."
    )
    parser.add_argument(
        "--name1",
        dest="inName1",
        required=True,
        help="Name for dataset 1."
    )
    parser.add_argument(
        "--prefix1",
        dest="inPrefix1",
        required=False,
        help="Prefix added to UJC for dataset 1 (after TranD prefix)."
    )
    parser.add_argument(
        "-2",
        "--dataset2",
        dest="inFile2",
        required=True,
        help="UJC key file for dataset 2."
    )
    parser.add_argument(
        "--name2",
        dest="inName2",
        required=True,
        help="Name for dataset 2."
    )
    parser.add_argument(
        "--prefix2",
        dest="inPrefix2",
        required=False,
        help="Prefix added to UJC for dataset 2 (after TranD prefix)."
    )
    parser.add_argument(
        "-u",
        "--union",
        dest="inUnion",
        required=True,
        help="UJC key file for union annotation."
    )

    # Output data
    parser.add_argument(
        "-o",
        "--output",
        dest="outFile",
        required=True,
        help="Output CSV key file for union annotation."
    )

    args = parser.parse_args()
    return args

def main():
    # Get TranD UJC key files

    firstConsolDf1 = pd.read_csv(args.inFile1, low_memory=False)
    firstConsolDf2 = pd.read_csv(args.inFile2, low_memory=False)
    secondConsolDf = pd.read_csv(args.inUnion, low_memory=False)
    # firstConsolDf1 = pd.read_csv("/Volumes/blue/mcintyre/share/transcript_distance/human_analysis/TranD_consol_gene_hg38_refseq_2_ensembl/transcript_id_2_consolidation_id.csv", low_memory=False)
    # firstConsolDf2 = pd.read_csv("/Volumes/blue/mcintyre/share/transcript_distance/human_analysis/TranD_consol_gene_hg38_ensembl/transcript_id_2_consolidation_id.csv", low_memory=False)
    # secondConsolDf = pd.read_csv("/Volumes/blue/mcintyre/share/transcript_distance/human_analysis/TranD_pairwise_ensembl_refseq_union/transcript_id_2_consolidation_id.csv", low_memory=False)

    # Check columns - include compatability with previous TranD output/utilities
    if "ujc_id" in firstConsolDf1.columns:
        colName1 = "ujc_id"
    elif "consolidation_transcript_id" in firstConsolDf1.columns:
        colName1 = "consolidation_transcript_id"
    else:
        print("!!! ERROR: UJC key file of dataset1 not properly formatted.")
        exit()
    if "ujc_id" in firstConsolDf2.columns:
        colName2 = "ujc_id"
    elif "consolidation_transcript_id" in firstConsolDf2.columns:
        colName2 = "consolidation_transcript_id"
    else:
        print("!!! ERROR: UJC key file of dataset2 not properly formatted.")
        exit()
    if "ujc_id" in secondConsolDf.columns:
        colNameU = "ujc_id"
    elif "consolidation_transcript_id" in secondConsolDf.columns:
        colNameU = "consolidation_transcript_id"
    else:
        print("!!! ERROR: UJC key file of union not properly formatted.")
        exit()

    colLst = ['gene_id', 'transcript_id', 'consolidation_transcript_id', 'ujc_id']


    firstConsolDf1 = firstConsolDf1[[c for c in colLst if c in firstConsolDf1.columns]]
    firstConsolDf2 = firstConsolDf2[[c for c in colLst if c in firstConsolDf2.columns]]
    secondConsolDf = secondConsolDf[[c for c in colLst if c in secondConsolDf.columns]]

    # Get names and any prefix used after consolidation
    name1 = args.inName1
    name2 = args.inName2
    prefix1 = args.inPrefix1
    prefix2 = args.inPrefix2
    # name1 = "RefSeq"
    # name2 = "Ensembl"
    # prefix1 = "rs"
    # prefix2 = "ens"

    if prefix1 is not None and prefix2 is not None:
        # Add prefixes to first consolidation transcript_ids
        firstConsolDf1[colName1] = prefix1 + "_" + firstConsolDf1[colName1]
        firstConsolDf2[colName2] = prefix2 + "_" + firstConsolDf2[colName2]

    # Change column names
    firstConsolDf1 = firstConsolDf1.rename(columns={
        colName1: name1 + "_uniq_jxn_id",
        "transcript_id": name1 + "_transcript_id",
        "gene_id": name1 + "_gene_id"
    })
    firstConsolDf2 = firstConsolDf2.rename(columns={
        colName2: name2 + "_uniq_jxn_id",
        "transcript_id": name2 + "_transcript_id",
        "gene_id": name2 + "_gene_id"
    })
    secondConsolDf = secondConsolDf.rename(columns={
        colNameU: "union_uniq_jxn_id",
        "transcript_id": "individual_uniq_jxn_id"

    })

    # Merge original transcript id values to the combined union annotation key files
    mergeConsol1 = pd.merge(
        firstConsolDf1,
        secondConsolDf,
        how="outer",
        left_on=name1 + "_uniq_jxn_id",
        right_on="individual_uniq_jxn_id",
        indicator="merge_check1"
    )
    # Check that all unique junction ids of first dataset are in the union annotation (no left_only)
    if mergeConsol1["merge_check1"].value_counts()["left_only"] > 0:
        print("WARNING: Unexpected first merge")
    mergeConsol2 = pd.merge(
        firstConsolDf2,
        mergeConsol1,
        how="outer",
        left_on=name2 + "_uniq_jxn_id",
        right_on="individual_uniq_jxn_id",
        indicator="merge_check2"
    )
    # Check that all unique junction ids of second dataset are in the union annotation (no left_only)
    if mergeConsol2["merge_check2"].value_counts()["left_only"] > 0:
        print("WARNING: Unexpected second merge")
    # Check that all union annotations have been used
    if mergeConsol2.groupby(["merge_check1", "merge_check2"])["individual_uniq_jxn_id"].count()[("right_only","right_only")] > 0:
        print("WARNING: Unexpected second merge")

    # Get list of dataset1 and dataset2 transcript_id values for each union junction chain
    unionDf = mergeConsol2.groupby(["gene_id", "union_uniq_jxn_id"]).agg({
        name1 + "_transcript_id": lambda x: "" if len(x)==0 else "|".join(x.dropna().astype(str)),
        name1 + "_uniq_jxn_id": lambda x: "" if len(x)==0 else "|".join(x.dropna().drop_duplicates().astype(str)),
        name2 + "_transcript_id": lambda x: "" if len(x)==0 else "|".join(x.dropna().astype(str)),
        name2 + "_uniq_jxn_id": lambda x: "" if len(x)==0 else "|".join(x.dropna().drop_duplicates().astype(str)),
    }).reset_index()

    # Remove prefixes from unique junction id values
    if prefix1 is not None and prefix2 is not None:
        unionDf[name1 + "_uniq_jxn_id"] = unionDf[name1 + "_uniq_jxn_id"].str[len(prefix1)+1:]
        unionDf[name2 + "_uniq_jxn_id"] = unionDf[name2 + "_uniq_jxn_id"].str[len(prefix2)+1:]

    # Output final file
    unionDf.to_csv(args.outFile, index=False)


if __name__ == '__main__':
    # Parse command line arguments
    global args
    args = getOptions()
    main()
