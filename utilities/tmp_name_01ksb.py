# -*- coding: utf-8 -*-
"""
Created on Thu Jun  6 13:17:01 2024

@author: k.bankole
"""

import argparse
import pandas as pd
import re


def getOptions():
    # Parse command line arguments
    # parser = argparse.ArgumentParser(description="")

    # # Input data
    # parser.add_argument("-", "--", dest="", required=True, help="")

    # # Output data
    # parser.add_argument("-", "--", dest="", required=True, help="")

    # args = parser.parse_args()
    return args


def main():

    erpFile = "Z://SHARE/McIntyre_Lab/sex_specific_splicing/compare_fiveSpecies_er_vs_data_gtf/mel_sexdet_er_vs_data_pattern_file_FBgn0004652.csv"
    flagFile = "Z://SHARE/McIntyre_Lab/sex_specific_splicing/compare_fiveSpecies_er_vs_data_gtf/mel_sexdet_er_vs_data_flag_file_FBgn0004652.csv"

    erpDf = pd.read_csv(erpFile, low_memory=False)
    flagDf = pd.read_csv(flagFile, low_memory=False)

    tmpDf = flagDf[['geneID', 'exonRegion']].copy().drop_duplicates()
    tmpDf = tmpDf[tmpDf['exonRegion'].str.contains('ER')]
    tmpDf = tmpDf.groupby('geneID').apply('count').reset_index()

    numERDf = tmpDf.copy()
    numERDf.columns = ['geneID', 'numER']

    workingDf = erpDf[['jxnHash', 'geneID', 'ERP']].copy()

    ergDf = erpDf.groupby(['geneID', 'ERP']).agg({
        'jxnHash': set
    }).reset_index()

    ergDf['ERG'] = ergDf['geneID'] + ":" + ergDf['ERP']

    ergDf = ergDf[['ERG', 'jxnHash', 'ERP']]

    ergDf['numJxnHash'] = ergDf['jxnHash'].apply(len)

    ergDf['numER'] = ergDf['ERP'].apply(len)

    # TODO: Reminder that transcripts with data only exons are being treated the same
    # dataOnlyERPDf = ergDf[ergDf['numER'] > 22]
    # noDataOnlyDf = ergDf.drop(dataOnlyERPDf.index)

    dfDf = ergDf.copy()
    dfDf['flagHasDataOnlyExon'] = dfDf['numER'].apply(
        lambda x: 1 if x > 22 else 0)

    # Pattern discernment!
    # 1. flag transcripts with all exon regions in the gene (may need to change depending on what happens to the data onlys)
    # Pattern: all 1s

    noDataOnlyDf['flag_hasAllER'] = noDataOnlyDf['ERP'].apply(
        lambda x: 1 if all(char == 1 for char in x) else 0)

    # 2. flag transcripts with an exon skip (one missing ER between two present ERs)
    # Pattern 101

    noDataOnlyDf['flag_exonSkip'] = noDataOnlyDf['ERP'].apply(
        lambda x: 1 if re.search('1[0]+1', x) is not None else 0)
    noDataOnlyDf['numExonSkip'] = noDataOnlyDf['ERP'].apply(
        lambda x: len(re.findall('1[0]+1', x)))


if __name__ == '__main__':
    # Parse command line arguments
    global args
    args = getOptions()
    main()
