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
    parser = argparse.ArgumentParser(
        description="Count num jxnHash per ERP and describe ERP through "
                    "flags. Takes in \"ERP.csv\" file as input."
    )

    # Input data
    parser.add_argument(
        "-p",
        "--pattern-file",
        dest="erpFile",
        required=True,
        help="Location of ER pattern file"
    )

    # parser.add_argument(
    #     "-c",
    #     "--count-file",
    #     dest="countFile",
    #     required=False,
    #     help="Location of counts per jxnHash output from id_ujc"
    # )

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

    parser.add_argument(
        "-s",
        "--sample-ID",
        dest="sampleID",
        required=False,
        help="Optional SampleID. Will create a sampleID column in output."
    )

    args = parser.parse_args()
    return args


def main():

    # erpFile = "Z://SHARE/McIntyre_Lab/sex_specific_splicing/compare_fiveSpecies_er_vs_data_gtf/FBgn0004652_test_ERP.csv"
    # outdir = 'C://Users/knife/Desktop/Code Dumping Ground/mcintyre'

    # erpFile = "~/mclab/SHARE/McIntyre_Lab/sex_specific_splicing/compare_fiveSpecies_er_vs_data_gtf/mel_sexdet_er_vs_data_pattern_file_FBgn0004652.csv"

    erpFile = "/nfshome/k.bankole/mnt/exasmb.rc.ufl.edu-blue/mcintyre/share/transcript_ortholog/aln_ujc_erp_output/fiveSpecies_2_dsan1_ujc_sexDetSubset_er_vs_dsan_F_2_dsan1_ujc_updGeneID_ERP.csv"
    prefix = "prefix"
    sampleID = "dsan_F"
    outdir = ""

    erpFile = args.erpFile
    outdir = args.outdir
    prefix = args.prefix
    sampleID = args.sampleID

    inERPDf = pd.read_csv(erpFile, low_memory=False)

    numUniqJxnHash = inERPDf['jxnHash'].nunique()
    print("There are {} unique jxnHash in the ERP file.".format(numUniqJxnHash))

    # Group into ERGs using ERPs and do some edits to the DF

    erpDf = inERPDf.copy()
    erpDf['flagDataOnlyExon'] = erpDf['numDataOnlyExon'].apply(
        lambda x: 1 if x >= 1 else 0)

    # TODO: data only exons? are patterns that are the same but with a different number of dataonly exons different?
    erpDf = erpDf.groupby(['geneID', 'ERP', 'flagDataOnlyExon']).agg({
        'jxnHash': set,
        'strand': set,
        'seqname': set
    }).reset_index()

    singleStrandERP = erpDf['strand'].apply(lambda x: len(x) == 1)

    if not singleStrandERP.all():
        print("There are transcripts belonging to more than one strand. Quitting.")
        quit()
    else:
        erpDf['strand'] = erpDf['strand'].apply(
            lambda x: list(x)[0])

    singleChrERP = erpDf['seqname'].apply(lambda x: len(x) == 1)

    if not singleChrERP.all():
        print("There are transcripts belonging to more than one strand. Quitting.")
        quit()
    else:
        erpDf['seqname'] = erpDf['seqname'].apply(
            lambda x: list(x)[0])

    erpDf = erpDf[['ERP', 'flagDataOnlyExon',
                   'jxnHash', 'geneID', 'strand', 'seqname']]

    erpDf['numJxnHash'] = erpDf['jxnHash'].apply(len)

    erpDf['numAnnotatedER'] = erpDf['ERP'].apply(lambda x: x.count('1'))

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

    patternSeekDf['patternSeek'] = patternSeekDf['ERP'].str.split(
        '_').str[1]

    # Pattern discernment!
    # 1. flag transcripts with all exon regions in the gene and no reference exon regions
    patternSeekDf['flagNoSkip'] = patternSeekDf['patternSeek'].apply(
        lambda x: 1 if all(char == '1' for char in x) else 0)

    patternSeekDf['flagNovel'] = patternSeekDf['patternSeek'].apply(
        lambda x: 1 if all(char == '0' for char in x) else 0)

    # 2. flag transcripts with an exon skip (one missing ER between two present ERs)
    patternSeekDf['flagERSkip'] = patternSeekDf.apply(
        lambda x: 1 if re.search('(?<=1)+0+(?=1)+', x['patternSeek']) is not None else 0, axis=1)

    # Do not uncomment.
    # patternSeekDf['numExonSkip'] = patternSeekDf.apply(
    #     lambda x: len(re.findall('(?<=1)+0+(?=1)+', x['ERP'])), axis=1)

    # 3. 5' and 3' fragment (compared to the gene)
    patternSeekDf['flag5pFragment'] = patternSeekDf.apply(
        lambda x: 1 if re.search(
            "^1+0+$", x['patternSeek']) is not None else 0, axis=1)

    # if x['strand'] == '+' else '^0+1+$',
    # x['ERP'][0:numERDct[x['geneID']]]) is not None else 0, axis=1)

    patternSeekDf['flag3pFragment'] = patternSeekDf.apply(
        lambda x: 1 if re.search(
            '^0+1+$', x['patternSeek']) is not None else 0, axis=1)
    # if x['strand'] == '+' else "^1+0+$",
    # x['ERP'][0:numERDct[x['geneID']]]) is not None else 0, axis=1)

    # 4. internal fragment
    patternSeekDf['flagIntrnlFrgmnt'] = patternSeekDf.apply(
        lambda x: 1 if re.search('^0+1+0+$', x['patternSeek']) is not None else 0, axis=1)

    # 5. first/last ER present
    patternSeekDf['flagFirstER'] = patternSeekDf.apply(
        lambda x: 1 if re.search(
            '^1', x['patternSeek']) is not None else 0, axis=1)
    # if x['strand'] == '+' else "1$",
    # x['ERP'][0:numERDct[x['geneID']]]) is not None else 0, axis=1)

    patternSeekDf['flagLastER'] = patternSeekDf.apply(
        lambda x: 1 if re.search(
            '1$', x['patternSeek']) is not None else 0, axis=1)
    # if x['strand'] == '+' else "^1",
    # x['ERP'][0:numERDct[x['geneID']]]) is not None else 0, axis=1)

    # note to self:
    # 5' frag = only has 5' exons (only has exons on the 3' side)
    # 3' frag = only has 3' exons
    # internal frag: only has internal (0s on either end)

    # Check that the number of unique starting jxnHash matches the sum of numJxnHash column
    outNumCheck = patternSeekDf['numJxnHash'].sum()

    print("There are {} unique jxnHash in the output file.".format(outNumCheck))

    if numUniqJxnHash != outNumCheck:
        raise Exception(
            "The number of jxnHashes in the input does not match the number of jxnHashes in the input.")

    outPrefix = "{}/{}".format(outdir, prefix)

    # GTF1_vs_GTF2_flagERP

    outFile = "{}_flagERP.csv".format(outPrefix)

    outDf = patternSeekDf[[
        'ERP', 'flagDataOnlyExon', 'geneID', 'seqname', 'strand', 'numJxnHash', 'numAnnotatedER',
        'flagNoSkip', 'flagNovel', 'flagERSkip', 'flag5pFragment',
        'flag3pFragment', 'flagIntrnlFrgmnt', 'flagFirstER',
        'flagLastER']].copy()

    if sampleID:
        outDf['sampleID'] = sampleID

    outDf.to_csv(outFile, index=False)


if __name__ == '__main__':
    # Parse command line arguments
    global args
    args = getOptions()
    main()
