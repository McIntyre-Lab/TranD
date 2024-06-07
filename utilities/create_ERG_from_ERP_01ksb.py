# -*- coding: utf-8 -*-
"""
Created on Thu Jun  6 13:17:01 2024

@author: k.bankole
"""

import argparse
import pandas as pd
import re
import time


def getOptions():
    # Parse command line arguments
    parser = argparse.ArgumentParser(description="")

    # Input data
    parser.add_argument(
        "-p",
        "--pattern-file",
        dest="erpFile",
        required=True,
        help="Location of ER pattern file"
    )

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

    # # Output data
    parser.add_argument(
        "-o",
        "--outdir",
        dest="outdir",
        required=True,
        help="Output directory"
    )

    parser.add_argument(
        "-n",
        "--fileName",
        dest="prefix",
        required=True,
        help="Prefix for output files. Required."
    )

    args = parser.parse_args()
    return args


def main():

    # erpFile = "Z://SHARE/McIntyre_Lab/sex_specific_splicing/compare_fiveSpecies_er_vs_data_gtf/mel_sexdet_er_vs_data_pattern_file_FBgn0004652.csv"
    # flagFile = "Z://SHARE/McIntyre_Lab/sex_specific_splicing/compare_fiveSpecies_er_vs_data_gtf/mel_sexdet_er_vs_data_flag_file_FBgn0004652.csv"
    # countFile = "//exasmb.rc.ufl.edu/blue/mcintyre/share/rmg_lmm_dros_data/mel2dmel6_jxnHash_cnts_sumTR_FBgn0004652.csv"
    # outdir = 'C://Users/knife/Desktop/Code Dumping Ground/mcintyre'
    # prefix = "test"

    erpFile = args.erpFile
    flagFile = args.flagFile
    countFile = args.countFile
    outdir = args.outdir
    prefix = args.prefix

    erpDf = pd.read_csv(erpFile, low_memory=False)
    flagDf = pd.read_csv(flagFile, low_memory=False)

    # Create a dictionary of each gene and its number of reference exon regions
    tmpDf = flagDf[['geneID', 'exonRegion']].copy().drop_duplicates()
    tmpDf = tmpDf[tmpDf['exonRegion'].str.contains('ER')]
    tmpDf = tmpDf.groupby('geneID').apply('count').reset_index()

    numERDct = dict(zip(tmpDf['geneID'], tmpDf['exonRegion']))

    # Group into ERGs using ERPs and do some edits to the DF
    ergDf = erpDf.groupby(['geneID', 'ERP']).agg({
        'jxnHash': set,
        'strand': set
    }).reset_index()

    singleStrandERG = ergDf['strand'].apply(lambda x: len(x) == 1)

    if not singleStrandERG.all():
        print("There are transcripts belonging to more than one strand. Quitting.")
        quit()
    else:
        ergDf['strand'] = ergDf['strand'].apply(
            lambda x: list(x)[0])

    ergDf['ERG'] = ergDf['geneID'] + ":" + ergDf['ERP']

    ergDf = ergDf[['ERG', 'jxnHash', 'ERP', 'geneID', 'strand']]

    ergDf['numJxnHash'] = ergDf['jxnHash'].apply(len)

    ergDf['numER'] = ergDf['ERP'].apply(lambda x: x.count('1'))

    ergDf['flagDataOnlyExon'] = ergDf.apply(
        lambda x: 1 if len(x['ERP']) > numERDct[x['geneID']] else 0, axis=1)

    patternSeekDf = ergDf.copy()

    # erpLst = [
    #     '1'*22,
    #     '1'*23,
    #     '0'*21+'1',
    #     '0'*22+'1',
    #     '0'*19+'1'*3,
    #     '0'*19+'1'*4,
    #     '0'*10 + '101' + '1'*9,
    #     '1' + '0'*21 + '1',
    #     '1110111111101111101100',
    #     '11101111111011111011001',
    #     '0000000000000000111111',
    #     '00000000000000001111111',
    #     '0000000000000000000001',
    #     '00000000000000000000011',
    #     '1111111110000000000000',
    #     '1000000000000000000000',
    #     '0000000011110000000000',
    #     '00000000111100000000001',
    #     '00000000100000000000001',
    #     '0000000010000000000000'
    # ]

    # geneLst = ['FBgn0004652'] * len(erpLst)
    # strandLst = ['-'] * len(erpLst)

    # patternSeekDf = pd.DataFrame({
    #     'geneID': geneLst,
    #     'ERP': erpLst,
    #     'strand': strandLst
    # })

    # patternSeekDf['flagHasDataOnlyExon'] = patternSeekDf.apply(
    #     lambda x: 1 if len(x['ERP']) > numERDct[x['geneID']] else 0, axis=1)

    # TODO: to be consistent for now, i'm just ignoring the last dataOnly 1 flag for discerning patterns
    # do any patterns with dataOnly exons get any of these flags? which ones? Should I just ignore the last 1?

    # Pattern discernment!
    # 1. flag transcripts with all exon regions in the gene
    patternSeekDf['flagNoSkip'] = patternSeekDf['ERP'].apply(
        lambda x: 1 if all(char == '1' for char in x) else 0)

    # 2. flag transcripts with an exon skip (one missing ER between two present ERs)
    patternSeekDf['flagExonSkip'] = patternSeekDf.apply(
        lambda x: 1 if re.search('(?<=1)+0+(?=1)+', x['ERP'][0:numERDct[x['geneID']]]) is not None else 0, axis=1)

    patternSeekDf['numExonSkip'] = patternSeekDf.apply(
        lambda x: len(re.findall('(?<=1)+0+(?=1)+', x['ERP'][0:numERDct[x['geneID']]])), axis=1)

    # 3. 5' and 3' fragment (compared to the gene)
    patternSeekDf['flag5pFragment'] = patternSeekDf.apply(
        lambda x: 1 if re.search(
            "^1+0+$" if x['strand'] == '+' else '^0+1+$',
            x['ERP'][0:numERDct[x['geneID']]]) is not None else 0, axis=1)

    patternSeekDf['flag3pFragment'] = patternSeekDf.apply(
        lambda x: 1 if re.search(
            '^0+1+$' if x['strand'] == '+' else "^1+0+$",
            x['ERP'][0:numERDct[x['geneID']]]) is not None else 0, axis=1)

    # 4. internal fragment
    patternSeekDf['flagInternalFragment'] = patternSeekDf.apply(
        lambda x: 1 if re.search('^0+1+0+$', x['ERP'][0:numERDct[x['geneID']]]) is not None else 0, axis=1)

    # 5. first/last ER present
    patternSeekDf['flagFirstERPresent'] = patternSeekDf.apply(
        lambda x: 1 if re.search(
            '^1' if x['strand'] == '+' else "1$",
            x['ERP'][0:numERDct[x['geneID']]]) is not None else 0, axis=1)

    patternSeekDf['flagLastERPresent'] = patternSeekDf.apply(
        lambda x: 1 if re.search(
            '1$' if x['strand'] == '+' else "^1",
            x['ERP'][0:numERDct[x['geneID']]]) is not None else 0, axis=1)

    # note to self:
    # 5' frag = only has 5' exons (only has exons on the 3' side)
    # 3' frag = only has 3' exons
    # internal frag: only has internal (0s on either end)

    outFile = "{}/{}_data_ERP_based_ERG.csv".format(outdir, prefix)
    patternSeekDf[[
        'ERG', 'geneID', 'strand', 'numJxnHash', 'numER', 'flagDataOnlyExon',
        'flagNoSkip', 'flagExonSkip', 'numExonSkip', 'flag5pFragment',
        'flag3pFragment', 'flagInternalFragment', 'flagFirstERPresent',
        'flagLastERPresent']].to_csv(outFile, index=False)

    if countFile:
        alphatic = time.perf_counter()

        countDf = pd.read_csv(countFile, low_memory=False)

        omegatoc = time.perf_counter()
        print(
            f"Count file read complete! Took {(omegatoc-alphatic):0.4f} seconds.")

        # Count number xscript per ERG using count file
        xscript2ERGDf = ergDf[['ERG', 'jxnHash',
                               'geneID', 'numER']].explode('jxnHash')

        uniqPatternHashSet = set(xscript2ERGDf['jxnHash'])
        uniqCountHashSet = set(countDf['jxnHash'])

        erpOnlyHshLst = list(uniqPatternHashSet - uniqCountHashSet)
        countOnlyHshLst = list(uniqCountHashSet - uniqPatternHashSet)

        hashesInBoth = list(uniqPatternHashSet.intersection(uniqCountHashSet))

        countDf = countDf[countDf['jxnHash'].isin(hashesInBoth)].copy()
        xscript2ERGDf = xscript2ERGDf[xscript2ERGDf['jxnHash'].isin(
            hashesInBoth)].copy()

        mergeCountAndERGDf = pd.merge(
            xscript2ERGDf, countDf, on='jxnHash', how='outer', indicator='merge_check')

        if erpOnlyHshLst:
            print(
                "POSSIBLE ERROR: THERE ARE JXNHASHES THAT APPEAR IN THE ERP OUTPUT BUT NOT IN THE COUNT FILE")

            pd.Series(erpOnlyHshLst).to_csv(
                "{}/list_{}_ERG_count_inFile_only_jxnHash.txt".format(outdir, prefix), index=False, header=False)

        if countOnlyHshLst:
            print(
                "CAUTION: THERE ARE JXNHASHES THAT APPEAR IN THE COUNT FILE BUT NOT IN THE ERP OUTPUT")

            pd.Series(countOnlyHshLst).to_csv(
                "{}/list_{}_ERG_count_countFile_only_jxnHash.txt".format(outdir, prefix), index=False, header=False)

        # if (mergeCountAndERGDf['merge_check'] == 'left_only').any():
        #     print(
        #         "POSSIBLE ERROR: THERE ARE JXNHASHES THAT APPEAR IN THE ERP OUTPUT BUT NOT IN THE COUNT FILE")
        #     # quit()
        # if (mergeCountAndERGDf['merge_check'] == 'right_only').any():
        #     print(
        #         "CAUTION: THERE ARE JXNHASHES THAT APPEAR IN THE COUNT FILE BUT NOT IN THE ERP OUTPUT")
        #     # quit()

        mergeCountAndERGDf = mergeCountAndERGDf.drop('merge_check', axis=1)
        ergCountDf = mergeCountAndERGDf.groupby(['sample', 'geneID', 'ERG']).agg({
            'jxnHash': set,
            'numTranscripts': sum,
            'numER': max
        }).reset_index()

        ergCountDf['numJxnHash'] = ergCountDf['jxnHash'].apply(len)

        # TODO: need new name
        outCountFile = "{}/{}_ERG_count.csv".format(outdir, prefix)
        ergCountDf[[
            'sample', 'geneID', 'ERG', 'numJxnHash', 'numTranscripts',
            'numER'
        ]].to_csv(outCountFile, index=False)


if __name__ == '__main__':
    # Parse command line arguments
    global args
    args = getOptions()
    main()
