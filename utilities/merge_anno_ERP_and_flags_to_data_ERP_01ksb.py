#!/usr/bin/env python

import argparse
import pandas as pd
import numpy as np


def getOptions():
    # Parse command line arguments
    parser = argparse.ArgumentParser(description="")

    # parser.add_argument(
    #     "-p",
    #     "--pattern-file",
    #     dest="erpFile",
    #     required=True,
    #     help="Location of ER pattern file"
    # )

    # parser.add_argument(
    #     "-c",
    #     "--count-file",
    #     dest="countFile",
    #     required=False,
    #     help="Location of counts per jxnHash"
    # )

    # # # Output data
    # parser.add_argument(
    #     "-o",
    #     "--outdir",
    #     dest="outdir",
    #     required=True,
    #     help="Output directory"
    # )

    # parser.add_argument(
    #     "-n",
    #     "--data-filename",
    #     dest="fileName",
    #     required=True,
    #     help="Name of data GTF for output files. Required."
    # )

    # parser.add_argument(
    #     "-x",
    #     "--prefix",
    #     dest="prefix",
    #     required=False,
    #     help="Prefix for output files."
    # )

    args = parser.parse_args()
    return args


def main():

    # inAnnoERPFile = "/nfshome/k.bankole/mnt/exasmb.rc.ufl.edu-blue/mcintyre/share/sex_specific_splicing/fiveSpecies_annotations/fiveSpecies_2_dmel6_ujc_er_vs_fiveSpecies_2_dmel6_ujc_ERP.csv"
    # inDataERPFile = "/nfshome/k.bankole/mnt/exasmb.rc.ufl.edu-blue/mcintyre/share/sex_specific_splicing/compare_fiveSpecies_er_vs_data_gtf/fiveSpecies_2_dmel6_ujc_er_vs_mel_2_dmel6_uniq_jxnHash_ERP.csv"
    # in5SpFlagFile = "/nfshome/k.bankole/mnt/exasmb.rc.ufl.edu-blue/mcintyre/share/sex_specific_splicing/fiveSpecies_annotations/flag_fiveSpecies_2_dmel6_ujc.csv"

    inAnnoERPFile = "//exasmb.rc.ufl.edu/blue/mcintyre/share/sex_specific_splicing/fiveSpecies_annotations/fiveSpecies_2_dmel6_ujc_er_vs_fiveSpecies_2_dmel6_ujc_ERP.csv"
    inDataERPFile = "//exasmb.rc.ufl.edu/blue/mcintyre/share/sex_specific_splicing/compare_fiveSpecies_er_vs_data_gtf/mel_sexdet_fiveSpecies_2_dmel6_ujc_er_sexDetSubset_vs_mel_2_dmel6_uniq_jxnHash_sexDetSubset_ERP.csv"
    in5SpFlagFile = "//exasmb.rc.ufl.edu/blue/mcintyre/share/sex_specific_splicing/fiveSpecies_annotations/flag_fiveSpecies_2_dmel6_ujc.csv"

    inAnnoERPDf = pd.read_csv(inAnnoERPFile, low_memory=False)
    inDataERPDf = pd.read_csv(inDataERPFile, low_memory=False)
    in5SpFlagDf = pd.read_csv(in5SpFlagFile, low_memory=False)

    uniqERPGeneSet = set(inAnnoERPDf['geneID'])
    uniqDataGeneSet = set(inDataERPDf['geneID'])
    uniqFlagGeneSet = set(in5SpFlagDf['geneID'])

    geneInAllSet = set.intersection(
        uniqERPGeneSet, uniqDataGeneSet, uniqFlagGeneSet)

    erpRemovedGnLst = list(uniqERPGeneSet - geneInAllSet)
    flagRemovedGnLst = list(uniqFlagGeneSet - geneInAllSet)
    dataRemovedGnLst = list(uniqDataGeneSet - geneInAllSet)

    erpDf = inAnnoERPDf[inAnnoERPDf['geneID'].isin(geneInAllSet)].copy()
    dataDf = inDataERPDf[inDataERPDf['geneID'].isin(geneInAllSet)].copy()
    flagDf = in5SpFlagDf[in5SpFlagDf['geneID'].isin(geneInAllSet)].copy()
    flagDf['jxnHash'] = flagDf[[col for col in flagDf.columns if 'jxnHash' in col]]

    numDataJxnHash = dataDf['jxnHash'].nunique()
    numFlagJxnHash = flagDf['jxnHash'].nunique()
    numERPJxnHash = erpDf['jxnHash'].nunique()

    erpDf = erpDf[['jxnHash', 'geneID', 'ERP', 'IRER', 'strand']]

    merge1Df = pd.merge(erpDf, flagDf, how='outer', on=[
                        'geneID', 'jxnHash'], indicator='merge_check')

    # TODO: error and success message
    if not (merge1Df['merge_check'] == 'both').all():
        print("Something")
        # quit()
    else:
        print('merge 1 success')
        merge1Df.drop('merge_check', axis=1, inplace=True)

    # TODO: this
    # merge1Df['jxnHash'].count()

    # dataDf = dataDf[['jxnHash', 'geneID', 'ERP']]

    merge2Df = pd.merge(dataDf, merge1Df, how='outer', on=[
                        'jxnHash'], indicator='merge_check', suffixes=['_data', '_anno'])

    quickOut = merge2Df[merge2Df['merge_check'] == 'both'].copy()

    quickOut = quickOut[['jxnHash', 'ERP_data', 'ERP_anno', 'geneID_data',
                         'geneID_anno', 'strand_data', 'strand_anno', 'IRER_data', 'IRER_anno']]

    quickOut.to_csv(
        'Z://SHARE/McIntyre_Lab/sex_specific_splicing/fiveSpecies_2_dmel6_test_compare_anno_vs_data_ERP.csv', index=False)


    # TODO: this
    # if refOnlyGnLst:
    #     pd.Series(refOnlyGnLst).to_csv(
    #         outPrefix + "list_{}_vs_{}_anno_only_genes.txt".format(erName, dataName), index=False, header=False)
    # if dataOnlyGnLst:
    #     pd.Series(dataOnlyGnLst).to_csv(
    #         outPrefix + "list_{}_vs_{}_data_only_genes.txt".format(erName, dataName), index=False, header=False)
if __name__ == '__main__':
    # Parse command line arguments
    global args
    args = getOptions()
    # main()
