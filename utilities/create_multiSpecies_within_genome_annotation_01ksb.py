#!/usr/bin/env python

import argparse
import glob
import os
import pandas as pd
import trand.io


def getOptions():
    # Parse command line arguments
    parser = argparse.ArgumentParser(description="Input genomeName as seen in mapping/ujc output, path to UJC GTF files all mapped to said genome, and"
                                     "outfile for the GTF (path and name). Outputs a file only contains unique jxnHashes "
                                     "from all species within one genome.")

    # Input data
    parser.add_argument("-gn",
                        "--genome-name",
                        dest="genomeName",
                        required=True,
                        help="Name of the coordinates that all the UJC indexes output are mapped to."
                        "Must match the file name (ex: example_2_{genomeName}_ujc.gtf")

    parser.add_argument("-ig",
                        "--input-gtfs",
                        dest="inFolder",
                        required=True,
                        help="Path to location of all UJC GTFs mapped to same genome.")

    # Output data
    parser.add_argument("-o",
                        "--output-file",
                        dest="outGTF",
                        required=True,
                        help="Path and name of GTF to output to (ex: "
                        "/path/to/output/fiveSpecies_2_dmel6_ujc.gtf")

    args = parser.parse_args()
    return args


def main():
    # Parse command line arguments

    # genomeName = "dmel6"
    # inFolder = "/nfshome/k.bankole/mnt/exasmb.rc.ufl.edu-blue/mcintyre/share/sex_specific_splicing/mapped_ujc_id_ujc_output"
    # outGTF = "/nfshome/k.bankole/mnt/exasmb.rc.ufl.edu-blue/mcintyre/share/sex_specific_splicing/fiveSpecies_annotations/fiveSpecies_2_dmel6_ujc_roz.gtf"

    # inFolder = r"\\exasmb.rc.ufl.edu\blue\mcintyre\share\sex_specific_splicing\mapped_ujc_id_ujc_output"

    genomeName = args.genomeName
    inFolder = args.inFolder
    outGTF = args.outGTF

    # Create dct {ujc gtf name:GTF DF}
    indexDfDct = dict()
    allHashSet = set()

    for file in glob.glob(inFolder + '/*_2_{}_roz_ujc.gtf'.format(genomeName)) + glob.glob(inFolder + '/*_2_{}_ujc_roz.gtf'.format(genomeName)):
        print(file)

        fileName = os.path.basename(file)

        annoName = fileName.split('2')[0].rstrip('_')
        inDf = trand.io.read_exon_data_from_file(file)

        print(len(allHashSet))

        indexDfDct[annoName] = inDf
        allHashSet.update(set(inDf['transcript_id'].unique().tolist()))
        print(len(allHashSet))

    print("Total Number of Unique UJCs across all annotations: {}".format(
        len(allHashSet)))

    test = set()
    for annoName, df in list(indexDfDct.items()):
        print(annoName)
        test.update(set(df['transcript_id'].unique().tolist()))

    len(set(test))

    annoName = list(indexDfDct.keys())[0]
    concatDf = indexDfDct[annoName].copy(deep=True)

    for annoName, df in list(indexDfDct.items())[1:]:
        concatDf = pd.concat(
            [concatDf, df]).drop_duplicates().reset_index(drop=True)

    print("Total Number of Unique UJCs in output: {}".format(
        concatDf['transcript_id'].nunique()))

    # Could use some extra verification but... seems like this is good?

    print("Writing GTF...")
    if os.path.isfile(outGTF):
        os.remove(outGTF)

    trand.io.write_gtf(data=concatDf, out_fhs={"gtf": outGTF}, fh_name="gtf")

    print("Complete!")


if __name__ == '__main__':
    global args
    args = getOptions()
    main()
