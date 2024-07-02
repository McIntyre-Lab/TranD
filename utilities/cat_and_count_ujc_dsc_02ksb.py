#!/usr/bin/env python

import argparse
import pandas as pd
import time


def getOptions():
    # Parse command line arguments
    parser = argparse.ArgumentParser(description="Reads in all samples and cats/counts the ujc DSC to contain only unique UJCs. "
                                     "so that it can go into another script to create a GTF with only unique UJCs.")

    # Input data
    parser.add_argument("-d",
                        "--design-file",
                        dest="desiFile",
                        required=True,
                        help="Design file with list of samples")

    parser.add_argument("-i",
                        "--input-directory",
                        dest="inDir",
                        required=True,
                        help="Location of all UJC dscrptn files")

    parser.add_argument("-gn",
                        "--genome-name",
                        dest="genomeName",
                        required=False,
                        help="Name of genome the samples were aligned to")

    parser.add_argument("-ff",
                        "--file-format",
                        dest="fileFormat",
                        required=True,
                        help="File format of the ujc dscrptn files. Describes where the sample ID is. "
                        "ex: sampleID_primary_ujc_dscrptn.csv. FileFormat = {}_primary_ujc_dscrptn.csv. "
                        "use a {} to replace where the sampleID would be.")

    parser.add_argument("-t",
                        "--tech-rep",
                        dest="techRep",
                        default="TechRep",
                        help="Name of the tech rep column in the design file. default: TechRep")
    # Output data
    parser.add_argument("-o",
                        "--output-file",
                        dest="outFile",
                        required=True,
                        help="Where to output unique UJC dscrptn files")

    args = parser.parse_args()
    return args


def main():
    print("Loading...", flush=True)
    alphatic = time.perf_counter()

    # desiFile = "/nfshome/k.bankole/mnt/exasmb.rc.ufl.edu-blue/mcintyre/share/rmg_lmm_dros_data/design_files/sample_bc_design_w_origDataPath_04amm.csv"
    # inDir = "/nfshome/k.bankole/mnt/exasmb.rc.ufl.edu-blue/mcintyre/share/rmg_lmm_dros_data/ujc_from_read_aln"
    # outFile = "/nfshome/k.bankole/Desktop/test_dsc/out_ujc_dsc.csv"
    # genomeName = "dmel6"
    # techRep = "TechRep"
    # fileFormat = "test"

    desiFile = args.desiFile
    inDir = args.inDir
    outFile = args.outFile
    genomeName = args.genomeName
    techRep = args.techRep
    fileFormat = args.fileFormat

    print("Reading design file...")
    dsnDf = pd.read_csv(desiFile)

    # Subsetting the design file for testing purposese only
    # dsnDf = dsnDf[dsnDf['sample'].apply(lambda x: 'mel' in x)]

    dscFileDct = {sampleID: "{}/{}.csv".format(
        inDir, fileFormat).format(sampleID) for sampleID in dsnDf['sampleID']}

    # dsnDf['sampleID'] = dsnDf['sample'] + "_TR" + dsnDf[techRep].astype(str)
    # dscFileDct = {sampleID:"{}/{}_2_{}_ujc_dscrptn.csv".format(inDir, sampleID, genomeName) for sampleID in dsnDf['sampleID']}

    toc = time.perf_counter()
    print(f"Complete! Took {toc-alphatic:0.4f} seconds.")
    tic = time.perf_counter()

    print("Total number of samples for input design file: ")
    print(len(dscFileDct))

    print()
    print("Reading desired files...", flush=True)
    dscDfDct = dict()
    allHashLst = []

    # for index, (sampleID, file) in enumerate(dscFileDct.items()):
    #     if index > 10:
    #         break
    for sampleID, file in dscFileDct.items():
        inDf = pd.read_csv(file, low_memory=False)

        inDf['sampleID'] = sampleID

        reorderDf = inDf[['jxnHash', 'jxnString', 'donorStart',
                          'acceptorEnd', 'chr', 'strand', 'sampleID']]

        print(sampleID + " counts:", flush=True)

        print("Number of total rows: {}".format(len(reorderDf['jxnHash'])))
        print("Number of unique jxnHash: {}".format(
            reorderDf['jxnHash'].nunique()))

        allHashLst = allHashLst + reorderDf['jxnHash'].unique().tolist()
        dscDfDct[sampleID] = reorderDf

    toc = time.perf_counter()
    print(f"Complete! Took {toc-tic:0.4f} seconds.", flush=True)
    tic = time.perf_counter()

    print("Total number of jxnHash across all samples:", len(allHashLst))
    print("Total number of unique jxnHash across all samples:", len(set(allHashLst)))
    print("Expected duplicate jxnHashes:",
          (len(allHashLst)-len(set(allHashLst))))

    print()
    print("Checking duplicated rows...", flush=True)
    # VERIFIED THAT NOT ALL DUPLICATE JXNHASHES WILL LEAD TO A FULL DUPLICATE ROW
    # (variation in the start and end of the transcript cause this). It's actually really cool.
    # There are some M/F transcripts from the same rep that only have a slight variation in their
    # 3' or 5' end.
    concatDf = pd.concat(dscDfDct.values())
    concatDf['dupe_row'] = concatDf.drop(columns='sampleID').duplicated()
    concatDf['dupe_jxnHash'] = concatDf.drop(
        columns=['sampleID'] + [col for col in concatDf.columns if 'dupe' in col]).duplicated(subset='jxnHash')

    toc = time.perf_counter()
    print(f"Complete! Took {toc-tic:0.4f} seconds.")
    tic = time.perf_counter()

    print("Number of duplicate jxnHashes in the concatenated dsc dataframe:",
          concatDf['dupe_jxnHash'].sum())
    print("Number of duplicate rows (ignoring sampleID column) in the concatenated dsc dataframe:",
          concatDf['dupe_row'].sum())
    print("Differences in the above two numbers are caused by UJCs with variation at the 3' or 5' end. aka: Duplicate rows = exact same transcript across multiple samples, Duplicate jxnHash = same UJC across multiple samples")

    print()
    print("Creating output file...", flush=True)
    # Group UJCs, keep earliest start and latest end
    groupDf = concatDf.groupby('jxnHash').agg({
        'jxnString': 'first',
        'donorStart': "min",
        'acceptorEnd': "max",
        'chr': 'first',
        'strand': 'first',
        'sampleID': 'unique'
    }).reset_index()

    groupDf['flag_multiSample'] = groupDf['sampleID'].apply(
        lambda x: len(x) > 1)
    groupDf.drop(columns='sampleID', inplace=True)

    # any(groupDf.duplicated())

    print("Number of unique jxnHash in output (duplicates removed):",
          groupDf['jxnHash'].nunique())
    groupDf.to_csv(outFile, index=False)

    omegatoc = time.perf_counter()
    print(f"Complete! Took {omegatoc-tic:0.4f} seconds.", flush=True)

    print()
    print(
        f"Entire operation took {omegatoc-alphatic:0.4f} seconds!", flush=True)


if __name__ == '__main__':
    # Parse command line arguments
    global args
    args = getOptions()
    main()
