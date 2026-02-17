#!/usr/bin/env python

import argparse
import glob
import os
import pandas as pd


def getOptions():
    # Parse command line arguments
    parser = argparse.ArgumentParser(description="Create flag file.")

    # Input data
    parser.add_argument("-a",
                        "--annotation-index",
                        dest="annoIndex",
                        required=True,
                        nargs="+",
                        help="Full path to annotation aligned to its own genome ujc xscript link.")

    parser.add_argument("-gn",
                        "--genome-name",
                        dest="genomeName",
                        required=True,
                        help="Name of the coordinates that all the UJC indexes output are mapped to."
                        "Must match the file name (ex: example_2_{genomeName}_ujc_xscript_link.csv.")

    parser.add_argument("-i",
                        "--input-directory",
                        dest="inDir",
                        required=True,
                        help="Path to location of all UJC indexes.")

    parser.add_argument("-gk",
                        "--gffcompare-gene-key",
                        dest="geneKey",
                        required=True,
                        help="Path to gffCompare gene key generated when creating the fiveSpecies annotation for genome"
                        )

    # Output data
    parser.add_argument("-o",
                        "--output-file",
                        dest="outFile",
                        required=True,
                        help="Path and name of file to output to (ex: "
                        "/path/to/output/genome_jxnHash_merge_flags.csv")

    args = parser.parse_args()
    return args


def main():
    # Parse command line arguments

    genomeName = "dmel6"
    inDir = "/nfshome/k.bankole/mnt/exasmb.rc.ufl.edu-blue/mcintyre/share/sex_specific_splicing/mapped_ujc_id_ujc_output"
    outFile = "/nfshome/k.bankole/mnt/exasmb.rc.ufl.edu-blue/mcintyre/share/sex_specific_splicing/test_fiveSpecies_merge/key_fiveSpecies_jxnHash_2_dmel6_jxnHash.csv"
    geneKey = "/nfshome/k.bankole/mnt/exasmb.rc.ufl.edu-blue/mcintyre/share/sex_specific_splicing/fiveSpecies_annotations/fiveSpecies_2_dmel6_ujc_gene_key.csv"
    annoIndex = ["/nfshome/k.bankole/mnt/exasmb.rc.ufl.edu-blue/mcintyre/share/sex_specific_splicing/annotation_id_ujc_output/dmel650_2_dmel6_ujc_xscript_link.csv"]

    genomeName = args.genomeName
    inDir = args.inDir
    outFile = args.outFile
    geneKey = args.geneKey
    annoIndex = args.annoIndex

    # print("Number of unique UJCs when fiveSpecies mapped to {}: {}".format(genomeName, len(geneDfr['{}_jxnHash'.format(genomeName)])))

    # Create dct {ujc gtf name:path to corresponding index}
    indexDfrDct = dict()
    allHashLst = []

    # Grab self-mapped file
    for file in annoIndex:
        fileName = os.path.basename(file)
        colName = fileName.split("_2_{}".format(genomeName))[
            0] + '_2_{}'.format(genomeName) + '_ujc'
        inDfr = pd.read_csv(file, low_memory=False)[['jxnHash']]

        inDfr['{}_jxnHash'.format(genomeName)] = inDfr['jxnHash']

        inDfr.columns = [colName + '_jxnHash', '{}_jxnHash'.format(genomeName)]

        print(colName)
        print("Num UJC transcripts: {}".format(
            len(inDfr['{}_jxnHash'.format(genomeName)])))
        print("Num UJCs when mapped to {}: {}".format(
            genomeName, inDfr['{}_jxnHash'.format(genomeName)].nunique()))
        print()

        inDfr = inDfr.groupby('{}_jxnHash'.format(genomeName)).agg(
            lambda x: list(set(x))).reset_index()

        allHashLst = allHashLst + \
            inDfr['{}_jxnHash'.format(genomeName)].unique().tolist()

        inDfr = inDfr[[colName + '_jxnHash', '{}_jxnHash'.format(genomeName)]]

        indexDfrDct[colName] = inDfr

    # Grab all files mapped to specified genome name
    for file in glob.glob(inDir + "/*_2_{}_noGeneID_ujc_xscript_link.csv".format(genomeName)):

        fileName = os.path.basename(file)
        colName = fileName.split("_2_{}".format(genomeName))[0]
        inDfr = pd.read_csv(file, low_memory=False)[
            ['transcriptID', 'jxnHash']]

        # fiveSpecies jxnHash = jxnHash on the  genomeName input (what it looks like in the annotation)
        inDfr.columns = [colName + '_jxnHash', '{}_jxnHash'.format(genomeName)]

        print(colName)
        print("Num UJC transcripts: {}".format(
            len(inDfr['{}_jxnHash'.format(genomeName)])))
        print("Num UJCs when mapped to {}: {}".format(
            genomeName, inDfr['{}_jxnHash'.format(genomeName)].nunique()))
        print()

        inDfr = inDfr.groupby('{}_jxnHash'.format(genomeName)).agg(
            lambda x: list(set(x))).reset_index()

        allHashLst = allHashLst + \
            inDfr['{}_jxnHash'.format(genomeName)].unique().tolist()

        inDfr = inDfr[[colName + '_jxnHash', '{}_jxnHash'.format(genomeName)]]

        indexDfrDct[colName] = inDfr

    totalRow = len(allHashLst)
    uniqHash = len(set(allHashLst))

    print("Num total UJC transcript rows: " + str(totalRow))
    print("Num total unique jxnHash when each annotation is mapped to {}: ".format(
        genomeName) + str(uniqHash))

    geneDfr = pd.read_csv(geneKey, low_memory=False)[
        ['transcript_id', 'output_gene_id']]
    geneDfr.columns = ['{}_jxnHash'.format(genomeName), 'geneID']

    # Start merge with first dataframe
    mergeDfr = geneDfr.copy(deep=True)
    # colName = list(indexDfrDct.keys())[0]
    # mergeDfrr['flag_{}'.format(colName)] = 1
    # print(mergeDfr['flag_{}'.format(colName)].sum())

    # Loop through all dataframes and merge on jxnHash
    for colName, dfr in list(indexDfrDct.items()):
        mergeDfr = pd.merge(mergeDfr, dfr, on=['{}_jxnHash'.format(
            genomeName)], how='outer', indicator='merge_check')

        if 'right_only' in mergeDfr['merge_check'].values:
            print("ERROR: Invalid Merge. Check your gene key.")
            exit()
        else:
            # Use merge_check to create flag for each UJC annotation
            mergeDfr['flag_{}'.format(colName)] = mergeDfr['merge_check'].apply(
                lambda x: 1 if x == 'right_only' or x == 'both' else 0)

            flagCol = [col for col in mergeDfr.columns if "flag" in col]
            mergeDfr[flagCol] = mergeDfr[flagCol].fillna(0)

            mergeDfr.drop(
                columns=['merge_check', colName + '_jxnHash'], inplace=True)

        # Verify flags (it works)
        # print(mergeDfr['flag_{}'.format(colName)].sum())

    firstCols = ['{}_jxnHash'.format(genomeName), 'geneID']
    otherCols = [col for col in mergeDfr.columns if col not in firstCols]
    outDfr = mergeDfr[firstCols + otherCols].copy(deep=True)

    # print()
    # print("Num output unique jxnHash: " + str(keyDfr['{}_jxnHash'.format(genomeName)].nunique()))

    # Check that the number of unique jxnHash in the index matches the
    # number of flagged jxnHash for each UJC annotation
    for colName, dfr in list(indexDfrDct.items()):
        numUniqHash = dfr['{}_jxnHash'.format(genomeName)].nunique()
        numFlagged = outDfr['flag_{}'.format(colName)].sum()

        print()
        print("Number " + colName + " unique jxnHash: " + str(numUniqHash))
        print("Number jxnHash flagged for " + colName + ": " + str(numFlagged))

        if (numUniqHash != numFlagged):
            print("ERROR: Number of unique hashes in input does not match"
                  "the total number of jxnHashes flagged for " + colName)
            exit()

    # Flag if jxnHash in all
    outDfr['numIn'] = outDfr.sum(axis=1, numeric_only=True)
    # outDfr['flag_all'] = outDfr['numIn'].apply(lambda x: 1 if x == 6 else 0)
    outDfr.drop(columns='numIn', inplace=True)

    # melBranchCol = [
    #     col for col in outDfr.columns if 'flag' in col and not 'dser' in col]
    # outDfr['numIn_melBranchOnly'] = outDfr[melBranchCol].sum(axis=1)
    # outDfr['flag_melBranchOnly'] = outDfr.apply(
    #     lambda x: 1 if x['numIn_melBranchOnly'] == 5 and x['flag_all'] != 1 else 0, axis=1)

    outDfr.drop(
        columns=[col for col in outDfr.columns if 'numIn' in col], inplace=True)
    outDfr.to_csv(outFile, index=False)

    # colLst = keyDfr.columns.tolist()
    # colLst.insert(0, colLst.pop(colLst.index('{}_jxnHash'.format(genomeName))))
    # keyDfr = keyDfr[colLst]

    # Output

    # Old Flag Code
    # Sort so jxnHash in all are at the top
    # outDfr = outDfr.sort_values(by='flag_all', ascending=False)
    # outDfr[outDfr.columns.difference(['jxnHash'])] = outDfr[outDfr.columns.difference(['jxnHash'])].astype(int)


if __name__ == '__main__':
    global args
    args = getOptions()
    main()
