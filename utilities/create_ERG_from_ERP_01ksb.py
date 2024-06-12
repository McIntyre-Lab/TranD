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
        "--data-filename",
        dest="fileName",
        required=True,
        help="Name of data GTF for output files. Required."
    )

    parser.add_argument(
        "-x",
        "--prefix",
        dest="prefix",
        required=False,
        help="Prefix for output files."
    )

    args = parser.parse_args()
    return args


def main():

    erpFile = "Z://SHARE/McIntyre_Lab/sex_specific_splicing/compare_fiveSpecies_er_vs_data_gtf/mel_sexdet_er_vs_data_pattern_file_FBgn0004652.csv"
    countFile = "//exasmb.rc.ufl.edu/blue/mcintyre/share/rmg_lmm_dros_data/mel2dmel6_jxnHash_cnts_sumTR_FBgn0004652.csv"
    outdir = 'C://Users/knife/Desktop/Code Dumping Ground/mcintyre'

    # erpFile = "~/mclab/SHARE/McIntyre_Lab/sex_specific_splicing/compare_fiveSpecies_er_vs_data_gtf/mel_sexdet_er_vs_data_pattern_file_FBgn0004652.csv"
    # countFile = "~/mnt/exasmb.rc.ufl.edu-blue/mcintyre/share/rmg_lmm_dros_data/mel2dmel6_jxnHash_cnts_sumTR_FBgn0004652.csv"

    prefix = "test"
    fileName = "gtf2"

    erpFile = args.erpFile
    countFile = args.countFile
    outdir = args.outdir
    fileName = args.fileName
    prefix = args.prefix

    inERPDf = pd.read_csv(erpFile, low_memory=False)

    # Group into ERGs using ERPs and do some edits to the DF
    erpDf = inERPDf.groupby(['geneID', 'ERP']).agg({
        'jxnHash': set,
        'strand': set
    }).reset_index()

    singleStrandERP = erpDf['strand'].apply(lambda x: len(x) == 1)

    if not singleStrandERP.all():
        print("There are transcripts belonging to more than one strand. Quitting.")
        quit()
    else:
        erpDf['strand'] = erpDf['strand'].apply(
            lambda x: list(x)[0])

    erpDf = erpDf[['ERP', 'jxnHash', 'geneID', 'strand']]

    erpDf['numJxnHash'] = erpDf['jxnHash'].apply(len)

    erpDf['numER'] = erpDf['ERP'].apply(lambda x: x.count('1'))

    patternSeekDf = erpDf.copy()

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

    # Pattern discernment!
    # 1. flag transcripts with all exon regions in the gene and no reference exon regions
    patternSeekDf['flag_noSkip'] = patternSeekDf['ERP'].apply(
        lambda x: 1 if all(char == '1' for char in x) else 0)

    patternSeekDf['flag_novel'] = patternSeekDf['ERP'].apply(
        lambda x: 1 if all(char == '0' for char in x) else 0)

    # 2. flag transcripts with an exon skip (one missing ER between two present ERs)
    patternSeekDf['flag_ERSkip'] = patternSeekDf.apply(
        lambda x: 1 if re.search('(?<=1)+0+(?=1)+', x['ERP']) is not None else 0, axis=1)

    # Do not uncomment.
    # patternSeekDf['numExonSkip'] = patternSeekDf.apply(
    #     lambda x: len(re.findall('(?<=1)+0+(?=1)+', x['ERP'])), axis=1)

    # 3. 5' and 3' fragment (compared to the gene)
    patternSeekDf['flag_5pFragment'] = patternSeekDf.apply(
        lambda x: 1 if re.search(
            "^1+0+$", x['ERP']) is not None else 0, axis=1)
    # if x['strand'] == '+' else '^0+1+$',
    # x['ERP'][0:numERDct[x['geneID']]]) is not None else 0, axis=1)

    patternSeekDf['flag_3pFragment'] = patternSeekDf.apply(
        lambda x: 1 if re.search(
            '^0+1+$', x['ERP']) is not None else 0, axis=1)
    # if x['strand'] == '+' else "^1+0+$",
    # x['ERP'][0:numERDct[x['geneID']]]) is not None else 0, axis=1)

    # 4. internal fragment
    patternSeekDf['flag_intrnlFrgmnt'] = patternSeekDf.apply(
        lambda x: 1 if re.search('^0+1+0+$', x['ERP']) is not None else 0, axis=1)

    # 5. first/last ER present
    patternSeekDf['flag_firstER'] = patternSeekDf.apply(
        lambda x: 1 if re.search(
            '^1', x['ERP']) is not None else 0, axis=1)
    # if x['strand'] == '+' else "1$",
    # x['ERP'][0:numERDct[x['geneID']]]) is not None else 0, axis=1)

    patternSeekDf['flag_lastER'] = patternSeekDf.apply(
        lambda x: 1 if re.search(
            '1$', x['ERP']) is not None else 0, axis=1)
    # if x['strand'] == '+' else "^1",
    # x['ERP'][0:numERDct[x['geneID']]]) is not None else 0, axis=1)

    # note to self:
    # 5' frag = only has 5' exons (only has exons on the 3' side)
    # 3' frag = only has 3' exons
    # internal frag: only has internal (0s on either end)

    # TODO: Check that the number of total unique jxnHash in pattern file matches
    # number of numJxnHash

    if prefix:
        outPrefix = "{}/{}_".format(outdir, prefix)
    else:
        outPrefix = "{}/".format(outdir)

    # GTF2_ERP_jxnHash_cnt

    outFile = outPrefix + "{}_ERP_jxnHash_cnt.csv".format(fileName)

    patternSeekDf[[
        'ERP', 'geneID', 'strand', 'numJxnHash', 'numER',
        'flag_noSkip', 'flag_novel', 'flag_ERSkip', 'flag_5pFragment',
        'flag_3pFragment', 'flag_intrnlFrgmnt', 'flag_firstER',
        'flag_lastER']].to_csv(outFile, index=False)

    # TODO: do this separately
    if countFile:
        alphatic = time.perf_counter()

        countDf = pd.read_csv(countFile, low_memory=False)

        omegatoc = time.perf_counter()
        print(
            f"Count file read complete! Took {(omegatoc-alphatic):0.4f} seconds.")

        # Count number xscript per ERG using count file
        xscript2ERPDf = erpDf[['ERP', 'jxnHash',
                               'geneID', 'numER', 'strand']].explode('jxnHash')

        uniqPatternHashSet = set(xscript2ERPDf['jxnHash'])
        uniqCountHashSet = set(countDf['jxnHash'])

        erpOnlyHshLst = list(uniqPatternHashSet - uniqCountHashSet)
        countOnlyHshLst = list(uniqCountHashSet - uniqPatternHashSet)

        hashesInBoth = list(uniqPatternHashSet.intersection(uniqCountHashSet))

        countDf = countDf[countDf['jxnHash'].isin(hashesInBoth)].copy()
        xscript2ERPDf = xscript2ERPDf[xscript2ERPDf['jxnHash'].isin(
            hashesInBoth)].copy()

        mergeCountAndERPDf = pd.merge(
            xscript2ERPDf, countDf, on='jxnHash', how='outer', indicator='merge_check')

        if erpOnlyHshLst:
            print(
                "POSSIBLE ERROR: THERE ARE JXNHASHES THAT APPEAR IN THE ERP OUTPUT BUT NOT IN THE COUNT FILE")

            pd.Series(erpOnlyHshLst).to_csv(
                outPrefix +
                "list_{}_ERP_jxnHash_cnt_hashes_in_ERPFile_only.txt".format(
                    fileName),
                index=False, header=False)

        if countOnlyHshLst:
            print(
                "CAUTION: THERE ARE JXNHASHES THAT APPEAR IN THE COUNT FILE BUT NOT IN THE ERP OUTPUT")

            outFile = outPrefix + \
                "list_{}_ERP_jxnHash_cnt_hashes_in_countFile_only.txt".format(
                    fileName)

            pd.Series(countOnlyHshLst).to_csv(
                outPrefix +
                "list_{}_ERP_jxnHash_cnt_hashes_in_countFile_only.txt".format(
                    fileName),
                index=False, header=False)

        # if (mergeCountAndERGDf['merge_check'] == 'left_only').any():
        #     print(
        #         "POSSIBLE ERROR: THERE ARE JXNHASHES THAT APPEAR IN THE ERP OUTPUT BUT NOT IN THE COUNT FILE")
        #     # quit()
        # if (mergeCountAndERGDf['merge_check'] == 'right_only').any():
        #     print(
        #         "CAUTION: THERE ARE JXNHASHES THAT APPEAR IN THE COUNT FILE BUT NOT IN THE ERP OUTPUT")
        #     # quit()

        mergeCountAndERPDf = mergeCountAndERPDf.drop('merge_check', axis=1)
        erpCountDf = mergeCountAndERPDf.groupby(['sample', 'geneID', 'ERG']).agg({
            'jxnHash': set,
            'numTranscripts': sum,
            'numER': max,
            'strand': set
        }).reset_index()

        singleStrandERPCnt = erpCountDf['strand'].apply(lambda x: len(x) == 1)

        if not singleStrandERPCnt.all():
            print("There are transcripts belonging to more than one strand. Quitting.")
            quit()
        else:
            erpCountDf['strand'] = erpCountDf['strand'].apply(
                lambda x: list(x)[0])

        erpCountDf['numJxnHash'] = erpCountDf['jxnHash'].apply(len)
        erpCountDf = erpCountDf.sort_values(by=['ERP', 'sample'])

        # TODO: need new name
        outCountFile = "{}/{}_ERP_count.csv".format(outdir, prefix)
        erpCountDf[[
            'sample', 'geneID', 'strand', 'ERP', 'numJxnHash', 'numTranscripts',
            'numER'
        ]].to_csv(outCountFile, index=False)


if __name__ == '__main__':
    # Parse command line arguments
    global args
    args = getOptions()
    main()
