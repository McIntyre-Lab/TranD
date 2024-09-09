#!/usr/bin/env python

import argparse
import pandas as pd
import time
import os
import warnings


def getOptions():
    # Parse command line arguments
    parser = argparse.ArgumentParser(
        description="Count num reads per ERP using read counts from "
                    "UJC output (ujc_count.csv) and \"ERP\" file. "
                    "Output will include counts per ERP in a stacked format. "
    )

    # Input data
    parser.add_argument(
        "-p",
        "--pattern-file",
        dest="inERPFile",
        required=True,
        help="ERP File"
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

    # inERPFile = "/nfshome/k.bankole/mnt/exasmb.rc.ufl.edu-blue/mcintyre/share/sex_specific_splicing/rlr_erp_output/fiveSpecies_2_dyak2_ujc_er_vs_dyak_data_2_dyak2_ujc_noMultiGene_ERP.csv"
    # inCntFile = "/nfshome/k.bankole/mnt/exasmb.rc.ufl.edu-blue/mcintyre/share/transcript_ortholog/dyak_data_2_dyak2_ujc_count.csv"

    inERPFile = "/nfshome/k.bankole/mnt/exasmb.rc.ufl.edu-blue/mcintyre/share/sex_specific_splicing/test_exon_segments_on_fru_dmel6/fiveSpecies_2_dmel6_ujc_Fru_er_vs_dmel_data_FBgn0004652_job_24_run_811_ujc_ERP.csv"
    inCntFile = "/nfshome/k.bankole/mnt/exasmb.rc.ufl.edu-blue/mcintyre/share/rmg_lmm_dros_data/ujc_byGene_output/dmel_data_FBgn0004652_job_24_run_811_ujc_count.csv"

    prefix = "fiveSpecies_2_dmel6_ujc_Fru_er_vs_dmel_data_FBgn0004652_job_24_run_811_ujc"
    outdir = "/nfshome/k.bankole/mnt/exasmb.rc.ufl.edu-blue/mcintyre/share/sex_specific_splicing/test_exon_segments_on_fru_dmel6"

    inERPFile = args.inERPFile
    inCntFile = args.inCntFile
    prefix = args.prefix
    outdir = args.outdir

    alphatic = time.perf_counter()

    inERPDf = pd.read_csv(inERPFile, low_memory=False)
    inCntDf = pd.read_csv(inCntFile, low_memory=False)

    uniqDataGeneSet = set(inCntDf['geneID'])
    uniqERPGeneSet = set(inERPDf['geneID'])

    # erpOnlyGnLst = list(uniqERPGeneSet - uniqDataGeneSet)
    dataOnlyGnLst = list(uniqDataGeneSet - uniqERPGeneSet)

    print("There should be", len(dataOnlyGnLst), "data only genes.")

    genesInBoth = list(uniqERPGeneSet.intersection(uniqDataGeneSet))

    erpDf = inERPDf[inERPDf['geneID'].isin(genesInBoth)].copy()
    cntDf = inCntDf[inCntDf['geneID'].isin(genesInBoth)].copy()

    if len(genesInBoth) != len(uniqDataGeneSet):
        warnings.warn("WARNING !!!! There are genes that are only in the count file. "
                      "If you did not subset the fiveSpecies ER annotation, "
                      "there is an issue. ")

    if len(genesInBoth) != len(uniqERPGeneSet):
        raise Exception("Error. There are genes in the count file that are not "
                        "in the ERP file. Be sure to subset the count file to "
                        "nonMultiGene jxnHash.")

    cntDf['numRead'] = cntDf['numRead'].astype(int)

    numHashPerSample = pd.DataFrame(cntDf.groupby(
        'sampleID')['jxnHash'].nunique().rename('startNum'))
    numReadPerSample = pd.DataFrame(cntDf.groupby(
        'sampleID')['numRead'].sum().rename('startNum'))

    toc = time.perf_counter()
    print(
        f"File read complete! Took {toc-alphatic:0.4f} seconds.")

    erpDf = erpDf[['jxnHash', 'ERP', 'geneID', 'strand', 'seqname']]

    # Merge ERP and counts
    erpMergeDf = pd.merge(erpDf, cntDf, on='jxnHash',
                          how='outer', indicator='merge_check')

    # Check that geneIDs are the same (I don't know why they wouldn't!)
    if all(erpMergeDf['geneID_x'] == erpMergeDf['geneID_y']):
        erpMergeDf['geneID'] = erpMergeDf['geneID_x']
        erpMergeDf.drop(['geneID_x', 'geneID_y'], inplace=True, axis=1)
    else:
        raise Exception("For some reason, a jxnHash has two different geneIDs in "
                        "the count file and the ERP file...")

    # Check that the merge did not add/remove jxnHash/reads per sample
    numHashPerSample['postERPMerge'] = erpMergeDf.groupby('sampleID').nunique()[
        'jxnHash']
    numReadPerSample['postERPMerge'] = erpMergeDf[[
        'numRead', 'sampleID']].groupby('sampleID').sum()

    if not (numHashPerSample['startNum'] == numHashPerSample['postERPMerge']).all():
        raise Exception("The merge of the ERP and Count Files led to a "
                        "different number of jxnHash per sample than before "
                        "the merge.")

    if not (numReadPerSample['startNum'] == numReadPerSample['postERPMerge']).all():
        raise Exception("The merge of the ERP and Count Files led to a "
                        "different number of reads per sample than before "
                        "the merge.")

    # Convert to unique on ERP and sum reads accross jxnHash
    erpCntDf = erpMergeDf.groupby(['sampleID', 'ERP', 'geneID']).agg({
        'jxnHash': 'size',
        'strand': set,
        'seqname': set,
        'numRead': sum
    }).reset_index()

    numHashPerSample['postERPGrp'] = erpMergeDf.groupby('sampleID').nunique()[
        'jxnHash']
    numReadPerSample['postERPGrp'] = erpMergeDf[[
        'numRead', 'sampleID']].groupby('sampleID').sum()

    if not (numHashPerSample['startNum'] == numHashPerSample['postERPGrp']).all():
        raise Exception("The merge of the ERP and Count Files led to a "
                        "different number of jxnHash per sample than before "
                        "the merge.")

    if not (numReadPerSample['startNum'] == numReadPerSample['postERPGrp']).all():
        raise Exception("The merge of the ERP and Count Files led to a "
                        "different number of reads per sample than before "
                        "the merge.")

    # Verify that total number of reads in converted Df matches counts in input
    if erpCntDf['numRead'].sum() != cntDf['numRead'].sum():
        raise Exception(
            "Error: Number of total reads in output does not match number of total reads in input.")

    # Make sure all ERPs are on a single strand and chr (dont know why they wouldnt be)
    singleStrandERP = erpCntDf['strand'].apply(lambda x: len(x) == 1)
    if not singleStrandERP.all():
        raise Exception(
            "There are ERPs belonging to more than one strand. Quitting.")
    else:
        erpCntDf['strand'] = erpCntDf['strand'].apply(
            lambda x: list(x)[0])

    singleChrERP = erpCntDf['seqname'].apply(lambda x: len(x) == 1)
    if not singleChrERP.all():
        raise Exception(
            "There are ERPs belonging to more than one seqname. Quitting.")
    else:
        erpCntDf['seqname'] = erpCntDf['seqname'].apply(
            lambda x: list(x)[0])

    erpCntDf = erpCntDf[['sampleID', 'ERP', 'geneID',
                         'strand', 'seqname', 'numRead']]

    print("ERP counts complete and verified!")

    outPrefix = f"{outdir}/{prefix}_"

    erpCntFile = outPrefix + "read_per_ERP.csv"

    erpCntDf.to_csv(
        erpCntFile, index=False)

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
