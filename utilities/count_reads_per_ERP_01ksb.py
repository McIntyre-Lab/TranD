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
                    "UJC output (ujc_count.csv) and \"infoERP\" file. "
                    "Output will include counts per ERP in a sample stacked format. "
    )

    # Input data
    parser.add_argument(
        "-i",
        "--info-file",
        dest="inInfoFile",
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

    inInfoFile = "/nfshome/k.bankole/mnt/exasmb.rc.ufl.edu-blue/mcintyre/share/sex_specific_splicing/rlr_erp_output/fiveSpecies_2_dyak2_ujc_er_vs_dyak_data_2_dyak2_ujc_noMultiGene_infoERP.csv"
    inCntFile = "/nfshome/k.bankole/mnt/exasmb.rc.ufl.edu-blue/mcintyre/share/transcript_ortholog/dyak_data_2_dyak2_ujc_count.csv"

    # inInfoFile = "/nfshome/k.bankole/mnt/exasmb.rc.ufl.edu-blue/mcintyre/share/sex_specific_splicing/rmg_erp_output/fiveSpecies_2_dmel6_ujc_er_vs_dmel_data_2_dmel6_ujc_noMultiGene_infoERP.csv"
    # inCntFile = "/nfshome/k.bankole/mnt/exasmb.rc.ufl.edu-blue/mcintyre/share/rmg_lmm_dros_data/dmel_data_2_dmel6_ujc_count.csv"

    prefix = "test"
    outdir = "/nfshome/k.bankole/Desktop/test_folder"
    # inERPFile = "/nfshome/k.bankole/mnt/exasmb.rc.ufl.edu-blue/mcintyre/share/sex_specific_splicing/test_exon_segments_on_fru_dmel6/fiveSpecies_2_dmel6_ujc_Fru_er_vs_dmel_data_FBgn0004652_job_24_run_811_ujc_ERP.csv"
    # inCntFile = "/nfshome/k.bankole/mnt/exasmb.rc.ufl.edu-blue/mcintyre/share/rmg_lmm_dros_data/ujc_byGene_output/dmel_data_FBgn0004652_job_24_run_811_ujc_count.csv"

    # prefix = "fiveSpecies_2_dmel6_ujc_Fru_er_vs_dmel_data_FBgn0004652_job_24_run_811_ujc"
    # outdir = "/nfshome/k.bankole/mnt/exasmb.rc.ufl.edu-blue/mcintyre/share/sex_specific_splicing/test_exon_segments_on_fru_dmel6"

    inInfoFile = args.inInfoFile
    inCntFile = args.inCntFile
    prefix = args.prefix
    outdir = args.outdir

    alphatic = time.perf_counter()

    inInfoDf = pd.read_csv(inInfoFile, low_memory=False)
    inCntDf = pd.read_csv(inCntFile, low_memory=False)

    uniqCntGeneSet = set(inCntDf['geneID'])
    uniqERPGeneSet = set(inInfoDf['geneID'])

    # erpOnlyGnLst = list(uniqERPGeneSet - uniqDataGeneSet)
    cntOnlyGnLst = list(uniqCntGeneSet - uniqERPGeneSet)

    print("There should be", len(cntOnlyGnLst), "data only genes.")

    genesInBoth = list(uniqERPGeneSet.intersection(uniqCntGeneSet))

    infoDf = inInfoDf[inInfoDf['geneID'].isin(genesInBoth)].copy()
    cntDf = inCntDf[inCntDf['geneID'].isin(genesInBoth)].copy()

    if len(genesInBoth) != len(uniqCntGeneSet):
        warnings.warn("WARNING !!!! There are genes that are only in the count file. "
                      "This could be due to TranD missing genes or the ER file being "
                      "a gene subset.")

    if len(genesInBoth) != len(uniqERPGeneSet):
        raise Exception("Error. There are genes in the count file that are not "
                        "in the infoERP file. Be sure to subset the count file to "
                        "noMultiGene jxnHash.")

    cntDf['numRead'] = cntDf['numRead'].astype(int)

    sumReadCnt = cntDf['numRead'].sum()
    cntGnDf = cntDf.groupby(['sampleID', 'geneID'])['numRead'].sum()

    infoDf = infoDf.fillna(0)
    infoDf['flagDataOnlyExon'] = inInfoDf['numDataOnlyExon'].apply(
        lambda x: 1 if x >= 1 else 0)

    numHashPerSample = pd.DataFrame(cntDf.groupby(
        'sampleID')['jxnHash'].nunique().rename('startNum'))
    numReadPerSample = pd.DataFrame(cntDf.groupby(
        'sampleID')['numRead'].sum().rename('startNum'))

    toc = time.perf_counter()
    print(
        f"File read complete! Took {toc-alphatic:0.4f} seconds.")

    infoDf = infoDf[['jxnHash', 'ERP', 'ERP_plus',
                     'geneID', 'flagDataOnlyExon', 'flagIR']]

    # Merge ERP and counts
    erpMergeDf = pd.merge(infoDf, cntDf, on='jxnHash',
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
    erpCntDf = erpMergeDf.groupby(['sampleID', 'geneID', 'ERP_plus']).agg({
        'jxnHash': 'size',
        'ERP': 'first',
        'flagDataOnlyExon': 'first',
        'flagIR': 'first',
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

    erpCntDf = erpCntDf[['sampleID', 'geneID', 'ERP', 'ERP_plus',
                         'flagDataOnlyExon', 'flagIR', 'numRead']]

    sumERPCnt = erpCntDf['numRead'].sum()
    erpGnDf = erpCntDf.groupby(['sampleID', 'geneID'])['numRead'].sum()

    if not (cntGnDf == erpGnDf).all() or sumERPCnt != sumReadCnt:
        raise Exception(
            "Error, the number of reads per gene/the number of total reads "
            "is not equivalent for input and output.")

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
