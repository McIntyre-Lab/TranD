#!/usr/bin/env python

import argparse
import pandas as pd
import time
import os
import warnings


def getOptions():
    # Parse command line arguments
    parser = argparse.ArgumentParser(
        description="Count num reads per ESP using read counts from "
                    "UJC output (ujc_count.csv) and \"ESP\" file. "
                    "Output will include counts per ESP in a stacked format. "
    )

    # Input data
    parser.add_argument(
        "-p",
        "--pattern-file",
        dest="inESPFile",
        required=True,
        help="ESP File"
    )

    parser.add_argument(
        "-c",
        "--count-file",
        dest="inCntFile",
        required=True,
        help="ujc_count File"
    )

    # Output data
    parser.add_argument(
        "-o",
        "--outdir",
        dest="outdir",
        required=True,
        help="Output directory. Must already exist."
    )

    parser.add_argument(
        "-x",
        "--prefix",
        dest="prefix",
        required=True,
        help="Prefix for output files. Required."
    )

    args = parser.parse_args()
    return args


def main():

    inESPFile = "/nfshome/k.bankole/mnt/exasmb.rc.ufl.edu-blue/mcintyre/share/sex_specific_splicing/test_exon_segments_on_fru_dmel6/fiveSpecies_2_dmel6_ujc_Fru_es_vs_dmel_data_FBgn0004652_job_24_run_811_ujc_ESP.csv"
    inCntFile = "/nfshome/k.bankole/mnt/exasmb.rc.ufl.edu-blue/mcintyre/share/rmg_lmm_dros_data/ujc_byGene_output/dmel_data_FBgn0004652_job_24_run_811_ujc_count.csv"

    prefix = "fiveSpecies_2_dmel6_ujc_Fru_es_vs_dmel_data_FBgn0004652_job_24_run_811_ujc"
    outdir = "/nfshome/k.bankole/mnt/exasmb.rc.ufl.edu-blue/mcintyre/share/sex_specific_splicing/test_exon_segments_on_fru_dmel6"

    inESPFile = args.inESPFile
    inCntFile = args.inCntFile
    prefix = args.prefix
    outdir = args.outdir

    alphatic = time.perf_counter()

    inESPDf = pd.read_csv(inESPFile, low_memory=False)
    inCntDf = pd.read_csv(inCntFile, low_memory=False)

    uniqDataGeneSet = set(inCntDf['geneID'])
    uniqESPGeneSet = set(inESPDf['geneID'])

    # espOnlyGnLst = list(uniqESPGeneSet - uniqDataGeneSet)
    dataOnlyGnLst = list(uniqDataGeneSet - uniqESPGeneSet)

    print("There should be", len(dataOnlyGnLst), "data only genes.")

    genesInBoth = list(uniqESPGeneSet.intersection(uniqDataGeneSet))

    espDf = inESPDf[inESPDf['geneID'].isin(genesInBoth)].copy()
    cntDf = inCntDf[inCntDf['geneID'].isin(genesInBoth)].copy()

    if len(genesInBoth) != len(uniqDataGeneSet):
        warnings.warn("WARNING !!!! There are genes that are only in the count file. "
                      "If you did not subset the fiveSpecies ER annotation, "
                      "there is an issue. ")

    if len(genesInBoth) != len(uniqESPGeneSet):
        raise Exception("Error. There are genes in the count file that are not "
                        "in the ESP file. Be sure to subset the count file to "
                        "nonMultiGene jxnHash.")

    cntDf['numRead'] = cntDf['numRead'].astype(int)

    numHashPerSample = pd.DataFrame(cntDf.groupby(
        'sampleID')['jxnHash'].nunique().rename('startNum'))
    numReadPerSample = pd.DataFrame(cntDf.groupby(
        'sampleID')['numRead'].sum().rename('startNum'))

    toc = time.perf_counter()
    print(
        f"File read complete! Took {toc-alphatic:0.4f} seconds.")

    espDf = espDf[['jxnHash', 'ESP', 'geneID', 'strand', 'seqname']]

    # Merge ESP and counts
    espMergeDf = pd.merge(espDf, cntDf, on='jxnHash',
                          how='outer', indicator='merge_check')

    # Check that geneIDs are the same (I don't know why they wouldn't!)
    if all(espMergeDf['geneID_x'] == espMergeDf['geneID_y']):
        espMergeDf['geneID'] = espMergeDf['geneID_x']
        espMergeDf.drop(['geneID_x', 'geneID_y'], inplace=True, axis=1)
    else:
        raise Exception("For some reason, a jxnHash has two different geneIDs in "
                        "the count file and the ESP file...")

    # Check that the merge did not add/remove jxnHash/reads per sample
    numHashPerSample['postESPMerge'] = espMergeDf.groupby('sampleID').nunique()[
        'jxnHash']
    numReadPerSample['postESPMerge'] = espMergeDf[[
        'numRead', 'sampleID']].groupby('sampleID').sum()

    if not (numHashPerSample['startNum'] == numHashPerSample['postESPMerge']).all():
        raise Exception("The merge of the ESP and Count Files led to a "
                        "different number of jxnHash per sample than before "
                        "the merge.")

    if not (numReadPerSample['startNum'] == numReadPerSample['postESPMerge']).all():
        raise Exception("The merge of the ESP and Count Files led to a "
                        "different number of reads per sample than before "
                        "the merge.")

    # Convert to unique on ESP and sum reads accross jxnHash
    espCntDf = espMergeDf.groupby(['sampleID', 'ESP', 'geneID']).agg({
        'jxnHash': 'size',
        'strand': set,
        'seqname': set,
        'numRead': sum
    }).reset_index()

    numHashPerSample['postESPGrp'] = espMergeDf.groupby('sampleID').nunique()[
        'jxnHash']
    numReadPerSample['postESPGrp'] = espMergeDf[[
        'numRead', 'sampleID']].groupby('sampleID').sum()

    if not (numHashPerSample['startNum'] == numHashPerSample['postESPGrp']).all():
        raise Exception("The merge of the ESP and Count Files led to a "
                        "different number of jxnHash per sample than before "
                        "the merge.")

    if not (numReadPerSample['startNum'] == numReadPerSample['postESPGrp']).all():
        raise Exception("The merge of the ESP and Count Files led to a "
                        "different number of reads per sample than before "
                        "the merge.")

    # Verify that total number of reads in converted Df matches counts in input
    if espCntDf['numRead'].sum() != cntDf['numRead'].sum():
        raise Exception(
            "Error: Number of total reads in output does not match number of total reads in input.")

    # Make sure all ESPs are on a single strand and chr (dont know why they wouldnt be)
    singleStrandESP = espCntDf['strand'].apply(lambda x: len(x) == 1)
    if not singleStrandESP.all():
        raise Exception(
            "There are ESPs belonging to more than one strand. Quitting.")
    else:
        espCntDf['strand'] = espCntDf['strand'].apply(
            lambda x: list(x)[0])

    singleChrESP = espCntDf['seqname'].apply(lambda x: len(x) == 1)
    if not singleChrESP.all():
        raise Exception(
            "There are ESPs belonging to more than one seqname. Quitting.")
    else:
        espCntDf['seqname'] = espCntDf['seqname'].apply(
            lambda x: list(x)[0])

    espCntDf = espCntDf[['sampleID', 'ESP', 'geneID',
                         'strand', 'seqname', 'numRead']]

    print("ESP counts complete and verified!")

    outPrefix = f"{outdir}/{prefix}_"

    espCntFile = outPrefix + "read_per_ESP.csv"

    espCntDf.to_csv(
        espCntFile, index=False)

    # if dataOnlyGnLst:
    #     pd.Series(dataOnlyGnLst).to_csv(
    #         outdir +
    #         "list_{}_2_{}_ujc_cnt_only_jxnHash.txt",
    #         index=False, header=False)

    omegatoc = time.perf_counter()

    print(f"Process Complete! Took {(omegatoc-alphatic):0.4f} seconds.")


if __name__ == '__main__':
    # Parse command line arguments
    global args
    args = getOptions()
    main()
