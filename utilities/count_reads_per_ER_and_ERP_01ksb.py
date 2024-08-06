#!/usr/bin/env python

import argparse
import pandas as pd
import time
import os
import warnings


def getOptions():
    # Parse command line arguments
    parser = argparse.ArgumentParser(
        description="Count num reads per ERP/ER using read counts from "
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

    # outdir = ""
    # prefix = None

    erpFileLst = args.erpFileLst
    erFileLst = args.erFileLst
    cntFileLst = args.cntFileLst
    genome = args.genome
    name = args.name
    sampleLst = args.sampleLst
    outdir = args.outdir
    prefix = args.prefix

    alphatic = time.perf_counter()

    inERFile = "/nfshome/k.bankole/mnt/exasmb.rc.ufl.edu-blue/mcintyre/share/transcript_ortholog/erp_output/fiveSpecies_2_dyak2_ujc_sexDetSubset_er_vs_dyak_data_2_dyak2_ujc_noMultiGene_flagER.csv"
    inERPFile = "/nfshome/k.bankole/mnt/exasmb.rc.ufl.edu-blue/mcintyre/share/transcript_ortholog/erp_output/fiveSpecies_2_dyak2_ujc_sexDetSubset_er_vs_dyak_data_2_dyak2_ujc_noMultiGene_ERP.csv"
    cntFile = "/nfshome/k.bankole/mnt/exasmb.rc.ufl.edu-blue/mcintyre/share/transcript_ortholog/ujc_species_gtf_combined_genes/dyak_data_2_dyak2_ujc_count_noMultiGene.csv"
    # genome = "dyak2"
    # name = "dyak_data"

    # Read in Files
    inERDf = pd.read_csv(inERFile, low_memory=False)
    inERPDf = pd.read_csv(inERPFile, low_memory=False)
    inCntDf = pd.read_csv(cntFile, low_memory=False)

    # verify that the ERP and flagER files have the same lis t of jxnHash
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

    erpDf = inERPDf[inERPDf['jxnHash'].isin(hashInBoth)].copy()
    erDf = inERDf[inERDf['jxnHash'].isin(hashInBoth)].copy()
    cntDf = inCntDf[inCntDf['jxnHash'].isin(hashInBoth)].copy()

    cntDf['numTranscripts'] = cntDf['numTranscripts'].astype(int)

    toc = time.perf_counter()
    print(
        f"File read complete! Took {toc-alphatic:0.4f} seconds.")

    ######### Part 1: Count Reads per ERP ########
    erpDf = erpDf[['jxnHash', 'ERP', 'geneID', 'strand', 'seqname']]

    # Merge ERP and counts
    erpMergeDf = pd.merge(erpDf, cntDf, on='jxnHash',
                          how='outer', indicator='merge_check')

    # Convert to unique on ERP and sum reads accross jxnHash
    erpCntDf = erpMergeDf.groupby(['source', 'ERP', 'geneID']).agg({
        'jxnHash': 'size',
        'strand': set,
        'seqname': set,
        'numTranscripts': sum
    }).reset_index()

    # Verify that total number of reads in converted Df matches counts in input
    if erpCntDf['numTranscripts'].sum() != cntDf['numTranscripts'].sum():
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

    erpCntDf = erpCntDf[['source', 'ERP', 'geneID',
                         'strand', 'seqname', 'numTranscripts']]
    erpCntDf = erpCntDf.rename(columns={'numTranscripts': 'numRead'})

    # pivot = pd.pivot_table(erpCntDf, index=['ERP','geneID','strand','seqname'], columns='source').reset_index()
    print("ERP counts complete and verified!")

    ######### Part 2: Count Reads per exon region (ER) ########

    erMergeDf = pd.merge(erDf, cntDf, on='jxnHash',
                         how='outer', indicator='merge_check')

    erMergeDf['numTranscripts'] = erMergeDf.apply(
        lambda x: x['numTranscripts'] if x['flagER'] == 1 else 0, axis=1)

    erCntDf = erMergeDf.groupby(['source', 'ER', 'geneID'], sort=False).agg({
        'flagER': sum,
        'numTranscripts': sum,
        'strand': set,
        'seqname': set,
        'lengthER': max
    }).reset_index()

    # TODO: need to figure out a verification step

    singleStrandER = erCntDf['strand'].apply(lambda x: len(x) == 1)
    if not singleStrandER.all():
        raise Exception(
            "There are ERs belonging to more than one strand. Quitting.")
    else:
        erCntDf['strand'] = erCntDf['strand'].apply(
            lambda x: list(x)[0])

    singleChrER = erCntDf['seqname'].apply(lambda x: len(x) == 1)
    if not singleChrER.all():
        raise Exception(
            "There are ERPs belonging to more than one chromosome. Quitting.")
    else:
        erCntDf['seqname'] = erCntDf['seqname'].apply(
            lambda x: list(x)[0])

    # Reorder columns to look nicer
    erCntDf = erCntDf[['source', 'ER', 'geneID',
                       'seqname', 'strand', 'numTranscripts']]

    print("ER counts complete and verified!")

    # TODO: output stuff
    outReadsPerERPDf = erpCntDf.copy()
    outReadsPerERDf = erCntDf.copy()

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
