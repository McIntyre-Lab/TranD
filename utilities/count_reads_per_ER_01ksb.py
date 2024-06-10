#!/usr/bin/env python

import argparse
import pandas as pd

def getOptions():
    # Parse command line arguments
    parser = argparse.ArgumentParser(description="")

    # # Input data
    # parser.add_argument("-", "--", dest="", required=True, help="")

    # # Output data
    # parser.add_argument("-", "--", dest="", required=True, help="")

    args = parser.parse_args()
    return args

def main():

    countFile = "~/mnt/exasmb.rc.ufl.edu-blue/mcintyre/share/rmg_lmm_dros_data/mel2dmel6_jxnHash_cnts_sumTR_FBgn0004652.csv"
    flagFile = "~/mclab/SHARE/McIntyre_Lab/sex_specific_splicing/compare_fiveSpecies_er_vs_data_gtf/mel_sexdet_er_vs_data_flag_file_FBgn0004652.csv"
    outdir  = ""
    prefix = ""
    
    inCountDf = pd.read_csv(countFile, low_memory=False)
    inFlagDf = pd.read_csv(flagFile, low_memory=False)
    
    uniqCountHashSet = set(inCountDf['jxnHash'])
    uniqFlagHashSet = set(inFlagDf['jxnHash'])

    countOnlyHshLst = list(uniqCountHashSet - uniqFlagHashSet)
    flagOnlyHashLst = list(uniqFlagHashSet - uniqCountHashSet)

    hashesInBoth = list(uniqFlagHashSet.intersection(uniqCountHashSet))

    countDf = inCountDf[inCountDf['jxnHash'].isin(hashesInBoth)].copy()
    flagDf = inFlagDf[inFlagDf['jxnHash'].isin(
        hashesInBoth)].copy()

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

    countDf['info'] = tuple(zip(countDf['sample'], countDf['numTranscripts']))    
    countDf = countDf[['jxnHash','info']].drop_duplicates()
    group = countDf.groupby('jxnHash').agg(list)
    
    mergeDf = pd.merge(group,inFlagDf,how='outer',on='jxnHash',indicator='merge_check')
    
    intmdDf = mergeDf.explode('info')
    intmdDf[['sample','numTranscripts']] = pd.DataFrame(intmdDf['info'].tolist(), index=intmdDf.index)

    intmdDf.drop(columns='info',inplace=True)

    
    # TODO: Check that the counts from the beginning match the end
    
    
if __name__ == '__main__':
    # Parse command line arguments
    global args
    # args = getOptions()
    main()
