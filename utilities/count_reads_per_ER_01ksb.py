#!/usr/bin/env python

import argparse
import pandas as pd


def getOptions():
    # Parse command line arguments
    parser = argparse.ArgumentParser(description="Placeholder")

    # INPUT
    parser.add_argument(
        "-f",
        "--flag-file",
        dest="flagFile",
        required=True,
        help="Location of data vs ref ER flag file"
    )

    parser.add_argument(
        "-c",
        "--count-file",
        dest="countFile",
        required=False,
        help="Location of counts per jxnHash"
    )

    # OUTPUT
    parser.add_argument(
        "-o",
        "--outdir",
        dest="outdir",
        required=True,
        help="Output directory"
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

    # countFile = "~/mnt/exasmb.rc.ufl.edu-blue/mcintyre/share/rmg_lmm_dros_data/mel2dmel6_jxnHash_cnts_sumTR_FBgn0004652.csv"
    # flagFile = "~/mclab/SHARE/McIntyre_Lab/sex_specific_splicing/compare_fiveSpecies_er_vs_data_gtf/mel_sexdet_er_vs_data_flag_file_FBgn0004652.csv"
    # outdir  = ""
    # prefix = ""

    # countFile = "//exasmb.rc.ufl.edu/blue/mcintyre/share/rmg_lmm_dros_data/mel2dmel6_jxnHash_cnts_sumTR_FBgn0004652.csv"
    # flagFile = "Z://SHARE/McIntyre_Lab/sex_specific_splicing/compare_fiveSpecies_er_vs_data_gtf/mel_sexdet_er_vs_data_flag_file_FBgn0004652.csv"
    # outdir = ""
    # prefix = ""

    flagFile = args.flagFile
    countFile = args.countFile
    outdir = args.outdir
    prefix = args.prefix

    inCountDf = pd.read_csv(countFile, low_memory=False)
    inFlagDf = pd.read_csv(flagFile, low_memory=False)

    uniqCountHashSet = set(inCountDf['jxnHash'])
    uniqFlagHashSet = set(inFlagDf['jxnHash'])

    countOnlyHshLst = list(uniqCountHashSet - uniqFlagHashSet)
    flagOnlyHashLst = list(uniqFlagHashSet - uniqCountHashSet)

    hashesInBoth = list(uniqFlagHashSet.intersection(uniqCountHashSet))

    if countOnlyHshLst:
        print(
            "POSSIBLE ERROR: THERE ARE JXNHASHES THAT APPEAR IN THE ERP OUTPUT BUT NOT IN THE COUNT FILE")

        pd.Series(countOnlyHshLst).to_csv(
            "{}/list_{}_ERG_count_inFile_only_jxnHash.txt".format(outdir, prefix), index=False, header=False)

    if flagOnlyHashLst:
        print(
            "CAUTION: THERE ARE JXNHASHES THAT APPEAR IN THE COUNT FILE BUT NOT IN THE ERP OUTPUT")

        pd.Series(flagOnlyHashLst).to_csv(
            "{}/list_{}_ERG_count_countFile_only_jxnHash.txt".format(outdir, prefix), index=False, header=False)

    countDf = inCountDf[inCountDf['jxnHash'].isin(hashesInBoth)].copy()
    flagDf = inFlagDf[inFlagDf['jxnHash'].isin(
        hashesInBoth)].copy()

    # z = countDf[countDf['jxnHash'] == '5a2f5b9accb2240be398a45b744c840ff2c6739013946781c2fc96950abe0caf'].copy()
    # y = flagDf[flagDf['jxnHash'] == "5a2f5b9accb2240be398a45b744c840ff2c6739013946781c2fc96950abe0caf"].copy()

    # z['info'] = tuple(zip(z['sample'], z['numTranscripts']))
    # y['info'] = tuple(zip(y['geneID'],y['strand'],y['exonRegion'],y['flag_ERPresent']))

    # zG = z.groupby('jxnHash').agg({'info':list}).reset_index()
    # yG = y.groupby('jxnHash').agg({'info':list}).reset_index()

    # inNumPerSample = countDf.groupby(['sample'])['numTranscripts'].sum()

    mergeDf = pd.merge(countDf, flagDf, on='jxnHash',
                       how='outer', indicator='merge_check')

    # TODO: check merge_check
    mergeDf.drop('merge_check', axis=1, inplace=True)

    mergeDf['actualNum'] = mergeDf['flag_ERPresent'] * \
        mergeDf['numTranscripts']

    outDf = mergeDf.groupby(['sample', 'exonRegion'], sort=False).agg({
        'geneID': set,
        'actualNum': sum,
        'strand': set,
        'flag_ERPresent': sum
    }).reset_index()

    singleGeneER = outDf['geneID'].apply(lambda x: len(x) == 1)
    if not singleGeneER.all():
        print("There are exon regions belonging to more than one gene. Quitting.")
        quit()
    else:
        outDf['geneID'] = outDf['geneID'].apply(
            lambda x: list(x)[0])

    singleStrandER = outDf['strand'].apply(lambda x: len(x) == 1)
    if not singleStrandER.all():
        print("There are exon regions belonging to more than one gene. Quitting.")
        quit()
    else:
        outDf['strand'] = outDf['strand'].apply(
            lambda x: list(x)[0])

    # outNumPerSample = outDf.groupby(['sample']).agg({
    #     'flag_ERPresent': 'sum'
    # })

    outDf = outDf[['geneID', 'strand', 'exonRegion', 'sample',
                   'actualNum']]

    outDf = outDf.rename({'actualNum': 'numTranscripts'}, axis=1)
    outFile = '{}/{}_data_ER_count.csv'.format(outdir, prefix)
    outDf[['geneID', 'strand', 'exonRegion', 'sample',
           'numTranscripts']].to_csv(outFile, index=False)

    # outDf = intmdDf.groupby(['sample', 'exonRegion'], sort=False).agg({
    #     'geneID': set,
    #     'actualNum': sum,
    #     'strand': set
    # }).reset_index()

    # singleGeneER = outDf['geneID'].apply(lambda x: len(x) == 1)

    # if not singleGeneER.all():
    #     print("There are exon regions belonging to more than one gene. Quitting.")
    #     quit()
    # else:
    #     outDf['geneID'] = outDf['geneID'].apply(
    #         lambda x: list(x)[0])

    # singleStrandER = outDf['strand'].apply(lambda x: len(x) == 1)

    # if not singleStrandER.all():
    #     print("There are exon regions belonging to more than one gene. Quitting.")
    #     quit()
    # else:
    #     outDf['strand'] = outDf['strand'].apply(
    #         lambda x: list(x)[0])

    # outDf = outDf.rename({'actualNum':'numTranscripts'}, axis=1)

    # outNumPerSample = outDf.groupby('sample')['numTranscripts'].sum()
    # TODO: Check that the counts from the beginning match the end
if __name__ == '__main__':
    # Parse command line arguments
    global args
    args = getOptions()
    main()
