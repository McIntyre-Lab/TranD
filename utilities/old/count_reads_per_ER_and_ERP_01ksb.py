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
                    "UJC output (ujc_count.csv) and \"ERP\"/\"flagER\" file. "
                    "Output will include counts per ERP/ER in a stacked format. "
    )

    # Input data
    parser.add_argument(
        "-p",
        "--pattern-file",
        dest="inERPFile",
        required=True,
        help="ERP File"
    )

    # parser.add_argument(
    #     "-f",
    #     "--flagER-file",
    #     dest="inERFile",
    #     required=True,
    #     help="flagER File"
    # )

    parser.add_argument(
        "-c",
        "--count-file",
        dest="inCntFile",
        required=True,
        help="ujc_count File"
    )

    parser.add_argument(
        "-g",
        "--genome",
        dest="genome",
        required=True,
        help="For the name of the output file (ex: name_2_genome_ERP_count.csv"
    )

    parser.add_argument(
        "-n",
        "--name",
        dest="name",
        required=True,
        help="Name for collective samples. For the name of the output file (ex: name_2_genome_ERP_count.csv)."
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

    # inERFile = "/nfshome/k.bankole/mnt/exasmb.rc.ufl.edu-blue/mcintyre/share/transcript_ortholog/erp_output/fiveSpecies_2_dyak2_ujc_sexDetSubset_er_vs_dyak_data_2_dyak2_ujc_noMultiGene_flagER.csv"
    # inERPFile = "/nfshome/k.bankole/mnt/exasmb.rc.ufl.edu-blue/mcintyre/share/transcript_ortholog/erp_output/fiveSpecies_2_dyak2_ujc_sexDetSubset_er_vs_dyak_data_2_dyak2_ujc_noMultiGene_ERP.csv"
    # inCntFile = "/nfshome/k.bankole/mnt/exasmb.rc.ufl.edu-blue/mcintyre/share/transcript_ortholog/ujc_species_gtf_combined_genes/dyak_data_2_dyak2_ujc_count_noMultiGene.csv"
    # genome = "dyak2"
    # name = "dyak_data"
    # outdir = "/nfshome/k.bankole/Desktop/test_folder"
    # prefix = None
    # prefix = "sexdet"

    # inERFile = args.inERFile
    inERPFile = args.inERPFile
    inCntFile = args.inCntFile
    genome = args.genome
    name = args.name
    outdir = args.outdir
    prefix = args.prefix

    alphatic = time.perf_counter()

    # Read in Files
    # inERDf = pd.read_csv(inERFile, low_memory=False)
    inERPDf = pd.read_csv(inERPFile, low_memory=False)
    inCntDf = pd.read_csv(inCntFile, low_memory=False)

    # verify that the ERP and flagER files have the same lis t of jxnHash
    # (if they don't there's an issue)
    # if len(inERPDf['jxnHash'].unique().tolist()) != inERDf['jxnHash'].nunique():
    #     raise Exception("Error. The jxnHashes in the ERP file are different from "
    #                     "the jxnHashes in the ER file. "
    #                     " Files: {} {}".format(inERPFile, inERFile))

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
    # erDf = inERDf[inERDf['jxnHash'].isin(hashInBoth)].copy()
    cntDf = inCntDf[inCntDf['jxnHash'].isin(hashInBoth)].copy()

    cntDf['numTranscripts'] = cntDf['numTranscripts'].astype(int)
    numHashPerSample = pd.DataFrame(cntDf.groupby('source').nunique()[
                                    'jxnHash'].rename('startNum'))
    numReadPerSample = pd.DataFrame(cntDf.groupby('source').sum()[
                                    'numTranscripts'].rename('startNum'))

    # sampleOvlpDf = cntDf.groupby('jxnHash').agg(list)[cntDf.groupby('jxnHash').agg(set)['source'].apply(lambda x: len(x) > 1)]

    toc = time.perf_counter()
    print(
        f"File read complete! Took {toc-alphatic:0.4f} seconds.")

    ######### Part 1: Count Reads per ERP ########
    erpDf = erpDf[['jxnHash', 'ERP', 'geneID', 'strand', 'seqname']]

    # Merge ERP and counts
    erpMergeDf = pd.merge(erpDf, cntDf, on='jxnHash',
                          how='outer', indicator='merge_check')

    # Check that the merge did not add/remove jxnHash/reads per sample
    numHashPerSample['postERPMerge'] = erpMergeDf.groupby('source').nunique()[
        'jxnHash']
    numReadPerSample['postERPMerge'] = erpMergeDf[[
        'numTranscripts', 'source']].groupby('source').sum()

    if not (numHashPerSample['startNum'] == numHashPerSample['postERPMerge']).all():
        raise Exception("The merge of the ERP and Count Files led to a "
                        "different number of jxnHash per sample than before "
                        "the merge.")

    if not (numReadPerSample['startNum'] == numReadPerSample['postERPMerge']).all():
        raise Exception("The merge of the ERP and Count Files led to a "
                        "different number of reads per sample than before "
                        "the merge.")

    # Convert to unique on ERP and sum reads accross jxnHash
    erpCntDf = erpMergeDf.groupby(['source', 'ERP', 'geneID']).agg({
        'jxnHash': 'size',
        'strand': set,
        'seqname': set,
        'numTranscripts': sum
    }).reset_index()

    # Check that the groupby did not add/remove jxnHash per sample
    numHashPerSample['postERPGrp'] = erpMergeDf.groupby('source').nunique()[
        'jxnHash']
    numReadPerSample['postERPGrp'] = erpMergeDf[[
        'numTranscripts', 'source']].groupby('source').sum()

    if not (numHashPerSample['startNum'] == numHashPerSample['postERPGrp']).all():
        raise Exception("The merge of the ERP and Count Files led to a "
                        "different number of jxnHash per sample than before "
                        "the merge.")

    if not (numReadPerSample['startNum'] == numReadPerSample['postERPGrp']).all():
        raise Exception("The merge of the ERP and Count Files led to a "
                        "different number of reads per sample than before "
                        "the merge.")

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
    print("ERP counts complete and verified!")

    # ######### Part 2: Count Reads per exon region (ER) ########

    # erMergeDf = pd.merge(erDf, cntDf, on='jxnHash',
    #                      how='outer', indicator='merge_check')

    # erMergeDf['numTranscripts'] = erMergeDf.apply(
    #     lambda x: x['numTranscripts'] if x['flagER'] == 1 else 0, axis=1)

    # # Check that the merge did not add/remove jxnHash per sample (reads are
    # # inherently counted multiple times due to there being more than 1 ER per read)
    # numHashPerSample['postERMerge'] = erMergeDf.groupby('source').nunique()[
    #     'jxnHash']

    # # numReadPerSample['postERMerge'] = erMergeDf[[
    # #     'numTranscripts', 'source']].groupby('source').sum()

    # if not (numHashPerSample['startNum'] == numHashPerSample['postERMerge']).all():
    #     raise Exception("The merge of the flagER and Count Files led to a "
    #                     "different number of jxnHash per sample than before "
    #                     "the merge.")

    # erCntDf = erMergeDf.groupby(['source', 'ER', 'geneID'], sort=False).agg({
    #     'flagER': sum,
    #     'numTranscripts': sum,
    #     'strand': set,
    #     'seqname': set,
    #     'lengthER': max
    # }).reset_index()

    # # TODO: need to figure out a verification step

    # singleStrandER = erCntDf['strand'].apply(lambda x: len(x) == 1)
    # if not singleStrandER.all():
    #     raise Exception(
    #         "There are ERs belonging to more than one strand. Quitting.")
    # else:
    #     erCntDf['strand'] = erCntDf['strand'].apply(
    #         lambda x: list(x)[0])

    # singleChrER = erCntDf['seqname'].apply(lambda x: len(x) == 1)
    # if not singleChrER.all():
    #     raise Exception(
    #         "There are ERPs belonging to more than one chromosome. Quitting.")
    # else:
    #     erCntDf['seqname'] = erCntDf['seqname'].apply(
    #         lambda x: list(x)[0])

    # # Reorder columns to look nicer
    # erCntDf = erCntDf[['source', 'ER', 'geneID',
    #                    'seqname', 'strand', 'numTranscripts']]

    # erCntDf = erCntDf.rename(columns={'numTranscripts': 'numRead'})

    # print("ER counts complete and verified!")

    # TODO: output stuff
    outReadsPerERPDf = erpCntDf.copy()
    # outReadsPerERDf = erCntDf.copy()

    if prefix:
        outPrefix = "{}/{}_".format(outdir, prefix)
    else:
        outPrefix = "{}/".format(outdir)

    readsPerERPFile = outPrefix + "{}_2_{}_ERP_cnts.csv".format(name, genome)

    outReadsPerERPDf.to_csv(
        readsPerERPFile, index=False)

    # readsPerERFile = outPrefix + "{}_2_{}_ER_cnts.csv".format(name, genome)
    # outReadsPerERDf.to_csv(
    #     readsPerERFile, index=False)

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

    # Code for converting from stacked to side-by-side + some stats stuff I was trying
    # erpCntWideDf = pd.pivot_table(erpCntDf, index=[
    #                               'ERP', 'geneID', 'strand', 'seqname'], columns='source', values=['numRead'], aggfunc=sum)
    # erpCntWideDf = erpCntWideDf.droplevel(0, axis=1).fillna(0)
    # erpCntWideDf.columns = [col + '_numRead' for col in erpCntWideDf.columns]

    # erpCntWideDf['sum_F_numRead'] = erpCntWideDf[[
    #     col for col in erpCntWideDf.columns if '_F_numRead' in col]].apply(sum, axis=1)

    # erpCntWideDf['sum_M_numRead'] = erpCntWideDf[[
    #     col for col in erpCntWideDf.columns if '_M_numRead' in col]].apply(sum, axis=1)

    # # statsTest = erpCntWideDf.copy()[['sum_F_numRead', 'sum_M_numRead']]

    # # statsTest['numDiff'] = statsTest['sum_F_numRead'] - statsTest['sum_M_numRead']
    # # statsTest = statsTest[statsTest['numDiff'] <= 574]

    # # avgDiff = statsTest['numDiff'].mean()
    # # stdDiff = statsTest['numDiff'].std()

    # # statsTest['z_score'] = (statsTest['numDiff'] - avgDiff) / stdDiff

    # erCntWideDf = pd.pivot_table(erCntDf, index=[
    #                              'ER', 'geneID', 'strand', 'seqname'], columns='source', values=['numRead'], aggfunc=sum)
    # erCntWideDf = erCntWideDf.droplevel(0, axis=1).fillna(0)
    # erCntWideDf.columns = [col + '_numRead' for col in erCntWideDf.columns]

    # erCntWideDf['sum_F_numRead'] = erCntWideDf[[
    #     col for col in erCntWideDf.columns if '_F_numRead' in col]].apply(sum, axis=1)

    # erCntWideDf['sum_M_numRead'] = erCntWideDf[[
    #     col for col in erCntWideDf.columns if '_M_numRead' in col]].apply(sum, axis=1)

    # # statsTest = erCntWideDf.copy()[['sum_F_numRead', 'sum_M_numRead']]

    # # statsTest['numDiff'] = statsTest['sum_F_numRead'] - statsTest['sum_M_numRead']
    # # statsTest = statsTest[statsTest['numDiff'] <= 544]

    # # avgDiff = statsTest['numDiff'].mean()
    # # stdDiff = statsTest['numDiff'].std()

    # # statsTest['z_score'] = (statsTest['numDiff'] - avgDiff) / stdDiff


if __name__ == '__main__':
    # Parse command line arguments
    global args
    args = getOptions()
    main()
