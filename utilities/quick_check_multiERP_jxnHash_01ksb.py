#!/usr/bin/env python

import argparse
import pandas as pd
import time
import os
import warnings


def getOptions():
    # Parse command line arguments
    parser = argparse.ArgumentParser(
        description="Count num reads per jxnHash using read counts from "
                    "UJC output. Combine across samples to create "
                    "side-by-side output. Needs gene_key created in GFFCompare "
                    "process to assign jxnHash to gene. For each "
                    "input file, input a list of all desired samples. "
                    "Output will include counts per ERP/ER split by sample, "
                    "where sample is in the column header."

    )

    # Input data
    parser.add_argument(
        "-c",
        "--count-file",
        dest="cntFileLst",
        required=True,
        nargs="+",
        help="One line list of counts per jxnHash output from id_ujc. Space separated."
    )

    parser.add_argument(
        "-gk",
        "--gene-key",
        dest="gnKeyLst",
        required=True,
        nargs="+",
        help="One line list of flagER files per sample. Space separated."
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
        "-n",
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

    # cntFileLst = [
    #     "/nfshome/k.bankole/mnt/exasmb.rc.ufl.edu-blue/mcintyre/share/transcript_ortholog/ujc_from_read_aln_samples/dyak_F_2_dyak2_ujc_count.csv",
    #     "/nfshome/k.bankole/mnt/exasmb.rc.ufl.edu-blue/mcintyre/share/transcript_ortholog/ujc_from_read_aln_samples/dyak_M_2_dyak2_ujc_count.csv"
    # ]

    # gnKeyLst = [
    #     "/nfshome/k.bankole/mnt/exasmb.rc.ufl.edu-blue/mcintyre/share/transcript_ortholog/gffcompare_read_aln_ujc/dyak_M_2_dyak2_ujc_gffcompare_gene_key.csv",
    #     "/nfshome/k.bankole/mnt/exasmb.rc.ufl.edu-blue/mcintyre/share/transcript_ortholog/gffcompare_read_aln_ujc/dyak_F_2_dyak2_ujc_gffcompare_gene_key.csv"
    # ]

    # sampleLst = ["dyak_F", "dyak_M"]
    # genome = "dyak2"
    # name = "dyak_data"

    # outdir = ""
    # prefix = None

    erpFileLst = [
        "//exasmb.rc.ufl.edu/blue/mcintyre/share/transcript_ortholog/aln_ujc_erp_output/fiveSpecies_2_dsan1_ujc_sexDetSubset_er_vs_dsan_F_2_dsan1_ujc_updGeneID_ERP.csv",
        "//exasmb.rc.ufl.edu/blue/mcintyre/share/transcript_ortholog/aln_ujc_erp_output/fiveSpecies_2_dsan1_ujc_sexDetSubset_er_vs_dsan_M_2_dsan1_ujc_updGeneID_ERP.csv"
    ]

    # gnKeyLst = [
    #     "//exasmb.rc.ufl.edu/blue/mcintyre/share/transcript_ortholog/gffcompare_read_aln_ujc/dsan_M_2_dsan1_ujc_gffcompare_gene_key.csv",
    #     "//exasmb.rc.ufl.edu/blue/mcintyre/share/transcript_ortholog/gffcompare_read_aln_ujc/dsan_F_2_dsan1_ujc_gffcompare_gene_key.csv"
    # ]

    sampleLst = ["dsan_F", "dsan_M"]
    genome = "dsan1"
    name = "dsan_data"

    outdir = ""
    prefix = None

    # gnKeyLst = args.gnKeyLst
    # erpFileLst = args.cntFileLst
    # sampleLst = args.sampleLst
    # genome = args.genome
    # name = args.name
    # outdir = args.outdir
    # prefix = args.prefix

    alphatic = time.perf_counter()

    # TODO: verification steps that I cannot think of rn

    # Loop through sample and get the three input files

    fileDfDct = dict()
    for sample in sampleLst:

        print(f"Reading {sample} files...")

        # Read in ERP file for sample

        # Read in flagER file for sample
        # inGnKey = [
        #     geneKey for geneKey in gnKeyLst
        #     if f"{sample}_2_{genome}" in geneKey
        #     and "gene_key.csv" in geneKey]

        # if len(inGnKey) == 1:
        #     inGnKey = inGnKey[0]
        # else:
        #     print(str(inGnKey))
        #     raise Exception(
        #         "There is either a duplicate input gene_key (two files with the same sample) "
        #         "or you have not input the correct file type for the -gk parameter.")

        # Read in ujc count file for sample
        inERPFile = [
            erpFile for erpFile in erpFileLst
            if f"{sample}_2_{genome}" in erpFile
            and "ERP.csv" in erpFile]

        if len(inERPFile) == 1:
            inERPFile = inERPFile[0]
        else:
            print(str(inERPFile))
            raise Exception(
                "There is either a duplicate input ujc count file (two files with the same sample) "
                "or you have not input the correct file type for the -c parameter.")

        # print("Gene Key:", inGnKey)
        print("ERP File:", inERPFile)

        # inGnKeyDf = pd.read_csv(inGnKey, low_memory=False)
        inERPDf = pd.read_csv(inERPFile, low_memory=False)

        # inGnKeyDf = inGnKeyDf[['transcript_id', 'output_gene_id']]
        # inGnKeyDf.columns = ['jxnHash', 'geneID']

        inERPDf = inERPDf[['jxnHash', 'ERP']]
        # verify that the ERP and flagER files have the same list of jxnHash
        # (if they don't there's an issue)
        # if len(inGnKeyDf['jxnHash'].unique().tolist()) != inCntDf['jxnHash'].nunique():
        #     raise Exception("Error. The jxnHashes in the gene key file are different from "
        #                     "the jxnHashes in the ujc count file. "
        #                     " Files: {} {}".format(inGnKey, inCntFile))

        # Only count for jxnHashes in both the count file and ERP file so that
        # script can be used on uneven subsets.
        # Output a really loud warning if any jxnHashes are removed
        # (counts will be off)

        # fileDfDct[sample] = [inGnKeyDf, inERPDf]
        fileDfDct[sample] = [inERPDf]

    toc = time.perf_counter()
    print(
        f"File read complete! Took {toc-alphatic:0.4f} seconds.")

    # For testing
    # sample = sampleLst[1]
    # test = fileDfDct[sample]
    # gnKeyDf = test[0]
    # cntDf = test[1]

    print("Merging gene into counts and stacking files...")
    stackedDf = pd.DataFrame()
    for sample, dfLst in fileDfDct.items():

        # gnKeyDf = dfLst[0]
        erpDf = dfLst[0]

        # mergeDf = pd.merge(gnKeyDf, erpDf, on=['jxnHash'],
        #                    how='outer', indicator='merge_check')
        erpDf['sampleID'] = sample

        # if (mergeDf['merge_check'] != 'both').all():
        #     raise Exception("I don't know how but when merging the gene key "
        #                     "and the count file, not all jxnHashes were the same. "
        #                     "I thought I already checked this...")
        # else:
        #     mergeDf.drop('merge_check', axis=1, inplace=True)

        if stackedDf.empty:
            stackedDf = erpDf
        else:
            stackedDf = pd.concat(
                [stackedDf, erpDf]).reset_index(drop=True)

    print("Creating flags...")
    numERPPerJxnHash = stackedDf.groupby('jxnHash')['ERP'].nunique()
    stackedDf['flagMultiERP'] = stackedDf['jxnHash'].map(
        numERPPerJxnHash > 1).astype(int)

    numSamplePerJxnHash = stackedDf.groupby(
        'jxnHash')['sampleID'].nunique()
    stackedDf['flagMultiSample'] = stackedDf['jxnHash'].map(
        numSamplePerJxnHash > 1).astype(int)

    multiSampleDf = stackedDf[stackedDf['flagMultiSample'] == 1]
    multiERPDf = multiSampleDf[multiSampleDf['flagMultiERP'] == 1]

    # Can be transitioned to wide mode if necessary (not complete)
    # pivot = geneAndCountDf.pivot(index=['geneID', 'jxnHash'],columns='sampleID')

    # # TODO: output stuff
    # outDf = geneAndCountDf[['sampleID', 'jxnHash', 'geneID',
    #                         'numRead', 'flagMultiGene', 'flagMultiSample']].copy()

    if prefix:
        outPrefix = "{}/{}_".format(outdir, prefix)
    else:
        outPrefix = "{}/".format(outdir)

    outFile = outPrefix + "{}_2_{}_jxnHash_cnts.csv".format(name, genome)

    outDf.to_csv(
        outFile, index=False)

    omegatoc = time.perf_counter()

    print(f"Process Complete! Took {(omegatoc-alphatic):0.4f} seconds.")


if __name__ == '__main__':
    # Parse command line arguments
    global args
    args = getOptions()
    main()

    # #san
    # sexDetLst = [
    #     "LOC120443850",
    #     "LOC120444305",
    #     "LOC120444548",
    #     "LOC120444712",
    #     "LOC120446026",
    #     "LOC120446047",
    #     "LOC120446856",
    #     "LOC120447068",
    #     "LOC120449354",
    #     "LOC120449495",
    #     "LOC120451005",
    #     "LOC120451624",
    #     "LOC120451673",
    #     "LOC120452969",
    #     "LOC120453140",
    #     "LOC120454172",
    #     "LOC120454598",
    #     "LOC120456630",
    #     "LOC120456785",
    #     "LOC120456786",
    #     "LOC120456871",
    #     "LOC120457037"
    # ]

    # #yak
    # sexDetLst = [
    #     "LOC6524166",
    #     "LOC6524656",
    #     "LOC6525459",
    #     "LOC6525973",
    #     "LOC6525974",
    #     "LOC6526726",
    #     "LOC6529542",
    #     "LOC6529909",
    #     "LOC6530293",
    #     "LOC6530404",
    #     "LOC6531347",
    #     "LOC6531813",
    #     "LOC6534435",
    #     "LOC6535120",
    #     "LOC6535575",
    #     "LOC6536551",
    #     "LOC6536981",
    #     "LOC6536985",
    #     "LOC6537555",
    #     "LOC6538122",
    #     "LOC6539897",
    #     "LOC6540023",
    # ]
