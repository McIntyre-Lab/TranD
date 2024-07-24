#!/usr/bin/env python

import argparse
import pandas as pd
import time
import os
import warnings


def getOptions():
    # Parse command line arguments
    parser = argparse.ArgumentParser(
        description="Count num reads per ERP/ERR using read counts from "
                    "UJC output and \"ERP\"/\"flagER\" file. For each "
                    "input file, input a list of all desired samples. "
                    "Output will include counts per ERP/ER split by sample, "
                    "where sample is in the column header."

    )

    # Input data
    parser.add_argument(
        "-p",
        "--pattern-file",
        dest="erpFileLst",
        required=True,
        nargs="+",
        help="One line list of ERP files per sample. Space separated."
    )

    parser.add_argument(
        "-f",
        "--flagER-file",
        dest="erFileLst",
        required=True,
        nargs="+",
        help="One line list of flagER files per sample. Space separated."
    )

    parser.add_argument(
        "-c",
        "--count-file",
        dest="cntFileLst",
        required=True,
        nargs="+",
        help="One line list of counts per jxnHash output from id_ujc. Space separated."
    )

    parser.add_argument(
        "-s",
        "--sampleID-list",
        dest="sampleLst",
        required=True,
        nargs="+",
        help="One line list of sample IDs. Space separated."
    )

    parser.add_argument(
        "-g",
        "--genome",
        dest="genome",
        required=True,
        help="Name of the coordinates that all the samples output are "
        "mapped to. Must match what is in the file name."
    )

    parser.add_argument(
        "-name",
        "--name",
        dest="name",
        required=True,
        help="Name for collective samples (ex: dmel_data)."
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
        required=False,
        help="Prefix for output files. Optional."
    )

    args = parser.parse_args()
    return args


def main():

    erpFileLst = [
        "/nfshome/k.bankole/mnt/exasmb.rc.ufl.edu-blue/mcintyre/share/transcript_ortholog/aln_ujc_erp_output/fiveSpecies_2_dsan1_ujc_sexDetSubset_er_vs_dsan_F_2_dsan1_ujc_updGeneID_ERP.csv",
        "/nfshome/k.bankole/mnt/exasmb.rc.ufl.edu-blue/mcintyre/share/transcript_ortholog/aln_ujc_erp_output/fiveSpecies_2_dsan1_ujc_sexDetSubset_er_vs_dsan_M_2_dsan1_ujc_updGeneID_ERP.csv"
    ]

    erFileLst = [
        "/nfshome/k.bankole/mnt/exasmb.rc.ufl.edu-blue/mcintyre/share/transcript_ortholog/aln_ujc_erp_output/fiveSpecies_2_dsan1_ujc_sexDetSubset_er_vs_dsan_F_2_dsan1_ujc_updGeneID_flagER.csv",
        "/nfshome/k.bankole/mnt/exasmb.rc.ufl.edu-blue/mcintyre/share/transcript_ortholog/aln_ujc_erp_output/fiveSpecies_2_dsan1_ujc_sexDetSubset_er_vs_dsan_M_2_dsan1_ujc_updGeneID_flagER.csv"
    ]

    cntFileLst = [
        "/nfshome/k.bankole/mnt/exasmb.rc.ufl.edu-blue/mcintyre/share/transcript_ortholog/ujc_from_read_aln_samples/dsan_F_2_dsan1_ujc_count.csv",
        "/nfshome/k.bankole/mnt/exasmb.rc.ufl.edu-blue/mcintyre/share/transcript_ortholog/ujc_from_read_aln_samples/dsan_M_2_dsan1_ujc_count.csv"
    ]

    sampleLst = ["dsan_F", "dsan_M"]
    genome = "dsan1"
    name = "dsan_data"
    outdir = ""
    prefix = None

    # erpFileLst = args.erpFileLst
    # flagERFileLst = args.flagERFileLst
    # cntFileLst = args.cntFileLst
    # genome = args.genome
    # name = args.name
    # sampleLst = args.sampleLst
    # outdir = args.outdir
    # prefix = args.prefix

    alphatic = time.perf_counter()

    # This is gonna suck.

    # TODO: verification steps that I cannot think of rn

    # Loop through sample and get the three input files

    print("Reading files...")
    fileDfDct = dict()
    for sample in sampleLst:
        # Read in ERP file for sample
        inERPFile = [
            erpFile for erpFile in erpFileLst
            if f"{sample}_2_{genome}" in erpFile
            and "ERP.csv" in erpFile]

        if len(inERPFile) == 1:
            inERPFile = inERPFile[0]
        else:
            raise Exception(
                "There is a duplicate input ERP file (two files with the same sample).")

        # Read in flagER file for sample
        inERFile = [
            erFile for erFile in erFileLst
            if f"{sample}_2_{genome}" in erFile
            and "flagER.csv" in erFile]

        if len(inERFile) == 1:
            inERFile = inERFile[0]
        else:
            raise Exception(
                "There is a duplicate input flagER file (two files with the same sample).")

        # Read in ujc count file for sample
        inCntFile = [
            cntFile for cntFile in cntFileLst
            if f"{sample}_2_{genome}" in cntFile
            and "ujc_count.csv" in cntFile]

        if len(inCntFile) == 1:
            inCntFile = inCntFile[0]
        else:
            raise Exception(
                "There is a duplicate input ujc count file (two files with the same sample).")

        print("Reading:", sample)
        print("ERP File:", inERPFile)
        print("ER File:", inERFile)
        print("Count File:", inCntFile)

        inERPDf = pd.read_csv(inERPFile, low_memory=False)
        inERDf = pd.read_csv(inERFile, low_memory=False)
        inCntDf = pd.read_csv(inCntFile, low_memory=False)

        # verify that the ERP and flagER files have the same list of jxnHash
        # (if they don't there's an issue)
        if len(inERPDf['jxnHash'].unique().tolist()) != inERDf['jxnHash'].nunique():
            raise Exception("Error. The jxnHashes in the ERP file are different from "
                            "the jxnHashes in the ER file. "
                            " Files: {} {}".format(inERPFile, inERFile))

        # Only count for jxnHashes in both the count file and ERP file so that
        # script can be used on uneven subsets.
        # Output a really loud warning if any jxnHashes are removed
        # (counts will be off)

        uniqERPHashSet = set(inERPDf['jxnHash'])
        uniqCntHashSet = set(inCntDf['jxnHash'])

        erpOnlyHshLst = list(uniqERPHashSet - uniqCntHashSet)
        cntOnlyHshLst = list(uniqCntHashSet - uniqERPHashSet)

        hashInBoth = list(uniqERPHashSet.intersection(uniqCntHashSet))

        if len(hashInBoth) != len(uniqCntHashSet):
            warnings.warn("WARNING !!!! There are jxnHashes that are only in the count file. "
                          "If you did not subset the fiveSpecies ER annotation, "
                          "there is an issue. ")

        if len(hashInBoth) != len(uniqERPHashSet):
            warnings.warn("WARNING !!!! There are jxnHashes that are only in the ERP file. "
                          "If you did not subset the UJC data files, "
                          "there is an issue. ")

        inERPDf = inERPDf[inERPDf['jxnHash'].isin(hashInBoth)].copy()
        inERDf = inERDf[inERDf['jxnHash'].isin(hashInBoth)].copy()
        inCntDf = inCntDf[inCntDf['jxnHash'].isin(hashInBoth)].copy()

        fileDfDct[sample] = [inERPDf, inERDf, inCntDf]

    toc = time.perf_counter()
    print(
        f"File read complete! Took {toc-alphatic:0.4f} seconds.")

    # # For testing
    # test = fileDfDct['dsan_F']
    # erpDf = test[0]
    # erDf = test[1]
    # cntDf = test[2]
    # sample = sampleLst[0]

    allSampleERPCntDf = pd.DataFrame()
    allSampleERCntDf = pd.DataFrame()
    for sample, dfLst in fileDfDct.items():

        # Get dataframes
        erpDf = dfLst[0]
        erDf = dfLst[1]
        cntDf = dfLst[2]
        numReadCol = sample + '_numRead'

        ######### Part 1: Count Reads per ERP ########
        erpDf = erpDf[['jxnHash', 'ERP', 'geneID', 'strand']]

        # Merge ERP and counts
        mergeDf = pd.merge(erpDf, cntDf, on='jxnHash',
                           how='outer', indicator='merge_check')
        mergeDf[numReadCol] = mergeDf['numTranscripts']

        # Convert to unique on ERP and sum reads accross jxnHash
        erpCntDf = mergeDf.groupby(['ERP', 'geneID']).agg({
            'jxnHash': 'size',
            'strand': set,
            numReadCol: sum
        }).reset_index()

        # Verify that total number of reads in converted Df matches counts in input
        if erpCntDf[numReadCol].sum() != cntDf['numTranscripts'].sum():
            raise Exception(
                "Error: Number of total reads in output does not match number of total reads in input.")

        # Make sure all ERPs are on a single strand (dont know why they wouldnt be)
        singleStrandERP = erpCntDf['strand'].apply(lambda x: len(x) == 1)
        if not singleStrandERP.all():
            raise Exception(
                "There are ERPs belonging to more than one strand. Quitting.")
        else:
            erpCntDf['strand'] = erpCntDf['strand'].apply(
                lambda x: list(x)[0])

        print(f"{sample} ERP counts complete and verified!")

        # Reorder columns to look nicer
        erpCntDf = erpCntDf[['ERP', 'geneID', 'strand', numReadCol]]

        ######### Part 2: Count Reads per exon region (ER) ########

        mergeDf = pd.merge(erDf, cntDf, on='jxnHash',
                           how='outer', indicator='merge_check')
        mergeDf[numReadCol] = mergeDf.apply(
            lambda x: x['numTranscripts'] if x['flagER'] == 1 else 0, axis=1)

        erCntDf = mergeDf.groupby(['ER', 'geneID'], sort=False).agg({
            'flagER': sum,
            numReadCol: sum,
            'strand': set,
        }).reset_index()

        # TODO: need to figure out a verification step

        singleStrandER = erCntDf['strand'].apply(lambda x: len(x) == 1)
        if not singleStrandER.all():
            raise Exception(
                "There are ERs belonging to more than one strand. Quitting.")
        else:
            erCntDf['strand'] = erCntDf['strand'].apply(
                lambda x: list(x)[0])

        print(f"{sample} ER counts complete and verified!")

        # Reorder columns to look nicer
        erCntDf = erCntDf[['ER', 'geneID', 'strand', numReadCol]]

        ######### Part 3: Merge samples together (keep samples in separate columns) ########

        if allSampleERPCntDf.empty:
            allSampleERPCntDf = erpCntDf
        else:
            allSampleERPCntDf = pd.merge(allSampleERPCntDf, erpCntDf, on=[
                                         'ERP', 'geneID', 'strand'], how='outer')

        readNumColLst = [
            col for col in allSampleERPCntDf.columns if 'numRead' in col]

        allSampleERPCntDf[readNumColLst] = \
            allSampleERPCntDf[readNumColLst].fillna(0)

        if allSampleERCntDf.empty:
            allSampleERCntDf = erCntDf
        else:
            allSampleERCntDf = pd.merge(allSampleERCntDf, erCntDf, on=[
                'ER', 'geneID', 'strand'], how='outer')

        readNumColLst = [
            col for col in allSampleERCntDf.columns if 'numRead' in col]

        allSampleERCntDf[readNumColLst] = \
            allSampleERCntDf[readNumColLst].fillna(0)

    # TODO: output stuff
    outReadsPerERPDf = allSampleERPCntDf.copy()
    outReadsPerERDf = allSampleERCntDf.copy()

    if prefix:
        outPrefix = "{}/{}_".format(outdir, prefix)
    else:
        outPrefix = "{}/".format(outdir)

    readsPerERPFile = outPrefix + "{}_2_{}_ERP_cnts.csv".format(name, genome)

    outReadsPerERPDf.to_csv(
        readsPerERPFile, index=False)

    readsPerERFile = outPrefix + "{}_2_{}_ER_cnts.csv".format(name, genome)
    outReadsPerERDf.to_csv(
        readsPerERFile, index=False)

    if erpOnlyHshLst:
        pd.Series(erpOnlyHshLst).to_csv(
            outPrefix +
            "list_{}_2_{}_ERP_only_jxnHash.txt".format(name, genome),
            index=False, header=False)

    if cntOnlyHshLst:
        pd.Series(cntOnlyHshLst).to_csv(
            outPrefix +
            "list_{}_2_{}_ujc_cnt_only_jxnHash.txt".format(name, genome),
            index=False, header=False)

    omegatoc = time.perf_counter()

    print(f"Process Complete! Took {(omegatoc-alphatic):0.4f} seconds.")


if __name__ == '__main__':
    # Parse command line arguments
    global args
    args = getOptions()
    main()
