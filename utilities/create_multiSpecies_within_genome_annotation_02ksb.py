#!/usr/bin/env python

import argparse
import glob
import os
import pandas as pd
import numpy as np
import trand.io


def getOptions():
    # Parse command line arguments
    parser = argparse.ArgumentParser(description="Input genomeName as seen in mapping/ujc output, "
                                     "path to UJC GTF files all mapped to said genome in one argument, and"
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
                        dest="inCat",
                        nargs="+",
                        required=True,
                        help="All paths to location of UJC GTFs mapped to same genome.")

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

    genomeName = "dmel6"
    inCatLst = ["/nfshome/k.bankole/mnt/exasmb.rc.ufl.edu-blue/mcintyre/share/sex_specific_splicing/annotation_id_ujc_output/dmel650_2_dmel6_ujc.gtf",
                "/nfshome/k.bankole/mnt/exasmb.rc.ufl.edu-blue/mcintyre/share/sex_specific_splicing/mapped_ujc_id_ujc_output/dsim202_2_dsim2_ujc_2_dmel6_noGeneID_ujc.gtf",
                "/nfshome/k.bankole/mnt/exasmb.rc.ufl.edu-blue/mcintyre/share/sex_specific_splicing/mapped_ujc_id_ujc_output/dsan11_2_dsan1_ujc_2_dmel6_noGeneID_ujc.gtf",
                "/nfshome/k.bankole/mnt/exasmb.rc.ufl.edu-blue/mcintyre/share/sex_specific_splicing/mapped_ujc_id_ujc_output/dsimWXD_2_dsim2_ujc_2_dmel6_noGeneID_ujc.gtf",
                "/nfshome/k.bankole/mnt/exasmb.rc.ufl.edu-blue/mcintyre/share/sex_specific_splicing/mapped_ujc_id_ujc_output/dyak21_2_dyak2_ujc_2_dmel6_noGeneID_ujc.gtf",
                "/nfshome/k.bankole/mnt/exasmb.rc.ufl.edu-blue/mcintyre/share/sex_specific_splicing/mapped_ujc_id_ujc_output/dser11_2_dser1_ujc_2_dmel6_noGeneID_ujc.gtf"]

    outGTF = "/nfshome/k.bankole/mnt/exasmb.rc.ufl.edu-blue/mcintyre/share/sex_specific_splicing/fiveSpecies_annotations/fiveSpecies_2_dmel6_ujc_roz.gtf"

    # genomeName = args.genomeName
    # inCatLst = args.inCat
    # outGTF = args.outGTF

    # TODO: keep track of number of exons for each transcript to check at the end

    gtfDct = dict()
    for file in inCatLst:

        fileName = os.path.basename(file)
        print("Reading file: {}".format(fileName))

        annoName = fileName.split('_2')[0].rstrip('_')
        df = trand.io.read_exon_data_from_file(file)
        gtfDct[annoName] = df

    rowLst = []

    gtfGrpDct = dict()
    for annoName, gtfDf in gtfDct.items():

        gtfDf['exonCoord'] = list(zip(gtfDf['start'], gtfDf['end']))

        gtfGrpDf = gtfDf.groupby(['seqname', 'strand', 'gene_id', 'transcript_id']).agg({
            'exonCoord': list
        }).reset_index()

        gtfGrpDct[annoName] = gtfGrpDf

    mergeDf = pd.DataFrame()
    suffix = 1
    for gtfGrpDf in list(gtfGrpDct.values()):

        if mergeDf.empty:
            mergeDf = gtfGrpDf
        else:
            mergeDf = pd.merge(mergeDf, gtfGrpDf, how='outer', on=[
                               'seqname', 'strand', 'transcript_id', 'gene_id'], indicator='merge_check', suffixes=[suffix, suffix+1])

            mergeDf['mergeCheck_' + str(suffix)] = mergeDf['merge_check']
            mergeDf = mergeDf.drop('merge_check', axis=1)

            suffix += 1

    rowLst = []

    # 14108
    for row in mergeDf.to_dict('records'):

        exonCoordLst = [value for key,
                        value in row.items() if 'exonCoord' in key]
        nonNanLst = [
            exonCoord for exonCoord in exonCoordLst if exonCoord is not np.nan]

        if len(nonNanLst) == 1:
            row['outExonLst'] = nonNanLst[0]

        else:

            # TODO: write this comment
            if any(len(exonCoord) != len(nonNanLst[0]) for exonCoord in nonNanLst):
                raise Exception("ruhroh")
            else:

                numExon = len(nonNanLst[0])
                exonLst = []
                # loop through the different set of coordinates that exist for each exon
                for exonCoordGrp in zip(*nonNanLst):

                    # if this exon has different coords across the species,
                    # find min start or last end and set exon to that.
                    if len(set(exonCoordGrp)) > 1:
                        exonPair = ()

                        if any(exonCoord[0] != exonCoordGrp[0][0] for exonCoord in exonCoordGrp):

                            exonLst.append(
                                (min(exonCoord[0] for exonCoord in exonCoordGrp), exonCoordGrp[0][1]))

                        elif any(exonCoord[1] != exonCoordGrp[0][1] for exonCoord in exonCoordGrp):

                            exonLst.append((exonCoordGrp[0][0], max(
                                exonCoord[1] for exonCoord in exonCoordGrp)))
                        else:
                            print("what happened???")

                    else:
                        exonLst.append(exonCoordGrp[0])

                row['outExonLst'] = exonLst

        rowLst.append(row)

    tmp = pd.DataFrame(rowLst)
    tmp = tmp.drop(columns=[col for col in tmp.columns if 'mergeCheck' in col])

    # kaboom = tmp[tmp['transcript_id'] ==
    #              "000ab5cb90829f133fd4e78f093a2d795bc47c83a7aaef93743560cf7c37f949"].copy()
    # kaboom = kaboom.dropna(axis=1)
    # kaboom = kaboom.explode(
    #     [col for col in kaboom.columns if 'exon' in col or 'Exon' in col])

    outDf = tmp.drop(
        columns=[col for col in tmp.columns if 'exonCoord' in col or 'mergeCheck' in col])

    outDf = outDf.explode('outExonLst')
    outDf[['start', 'end']] = pd.DataFrame(
        outDf['outExonLst'].to_list(), index=outDf.index)

    outTest = pd.DataFrame()
    outTest[['start', 'end']] = outDf['outExonLst'].apply(pd.Series)

    outDf = outDf[['seqname', 'strand', 'start',
                   'end', 'transcript_id', 'gene_id']]

    # verify exon counts.

    testDf = masterDf[masterDf['transcript_id'] ==
                      "000ab5cb90829f133fd4e78f093a2d795bc47c83a7aaef93743560cf7c37f949"]
    test2Df = gtfDf[gtfDf['transcript_id'] ==
                    "000ab5cb90829f133fd4e78f093a2d795bc47c83a7aaef93743560cf7c37f949"]

    print("Writing GTF...")
    if os.path.isfile(outGTF):
        os.remove(outGTF)

    trand.io.write_gtf(data=masterDf, out_fhs={
                       "gtf": "/nfshome/k.bankole/mclab/SHARE/McIntyre_Lab/sex_specific_splicing/one_gene_her_trand_test/test.gtf"}, fh_name="gtf")

    print("Complete!")


if __name__ == '__main__':
    global args
    args = getOptions()
    main()
