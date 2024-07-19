#!/usr/bin/env python

import argparse
import glob
import os
import pandas as pd
import numpy as np
import trand.io


def getOptions():
    # Parse command line arguments
    parser = argparse.ArgumentParser(description="Input path to UJC GTF files all mapped to "
                                     "one genome in one argument and "
                                     "outfile for the GTF (path and name). "
                                     "Outputs a file only contains unique jxnHashes "
                                     "from all species within one genome. If "
                                     "the same jxnHash exists across multiple "
                                     "annotations, the earliest start for the "
                                     "first exon and the latest end for the "
                                     "last exon will be chosen.")

    # Input data
    # parser.add_argument("-gn",
    #                     "--genome-name",
    #                     dest="genomeName",
    #                     required=True,
    #                     help="Name of the coordinates that all the UJC indexes output are mapped to."
    #                     "Must match the file name (ex: example_2_{genomeName}_ujc.gtf")

    parser.add_argument("-ig",
                        "--input-gtfs",
                        dest="inGTF",
                        nargs="+",
                        required=True,
                        help="All paths to location of UJC GTFs mapped to same genome in a "
                        "space separated format (all in one line).")

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
    # inGTFLst = ["/nfshome/k.bankole/mnt/exasmb.rc.ufl.edu-blue/mcintyre/share/sex_specific_splicing/annotation_id_ujc_output/dmel650_2_dmel6_ujc.gtf",
    #             "/nfshome/k.bankole/mnt/exasmb.rc.ufl.edu-blue/mcintyre/share/sex_specific_splicing/mapped_ujc_id_ujc_output/dsim202_2_dsim2_ujc_2_dmel6_noGeneID_ujc.gtf",
    #             "/nfshome/k.bankole/mnt/exasmb.rc.ufl.edu-blue/mcintyre/share/sex_specific_splicing/mapped_ujc_id_ujc_output/dsan11_2_dsan1_ujc_2_dmel6_noGeneID_ujc.gtf",
    #             "/nfshome/k.bankole/mnt/exasmb.rc.ufl.edu-blue/mcintyre/share/sex_specific_splicing/mapped_ujc_id_ujc_output/dsimWXD_2_dsim2_ujc_2_dmel6_noGeneID_ujc.gtf",
    #             "/nfshome/k.bankole/mnt/exasmb.rc.ufl.edu-blue/mcintyre/share/sex_specific_splicing/mapped_ujc_id_ujc_output/dyak21_2_dyak2_ujc_2_dmel6_noGeneID_ujc.gtf",
    #             "/nfshome/k.bankole/mnt/exasmb.rc.ufl.edu-blue/mcintyre/share/sex_specific_splicing/mapped_ujc_id_ujc_output/dser11_2_dser1_ujc_2_dmel6_noGeneID_ujc.gtf"]

    # outGTF = "/nfshome/k.bankole/mnt/exasmb.rc.ufl.edu-blue/mcintyre/share/sex_specific_splicing/fiveSpecies_annotations/fiveSpecies_2_dmel6_ujc_roz.gtf"

    inGTFLst = args.inGTF
    outGTF = args.outGTF

    allHashSet = set()
    # Read in all GTFs and store in dictionary
    gtfDct = dict()

    for gtf in inGTFLst:

        fileName = os.path.basename(gtf)
        print("Reading file: {}".format(fileName))

        annoName = fileName.split('_2')[0].rstrip('_')
        data = trand.io.read_exon_data_from_file(gtf)
        gtfDct[annoName] = data

        allHashSet.update(set(data['transcript_id'].unique().tolist()))

    print()
    print("GTF Read Complete!")
    print()

    print("Total Number of Unique UJCs across all annotations: {}".format(
        len(allHashSet)))

    # Group all GTF by anything unique to a transcript and store the exon coordinates in a list
    grpDct = dict()
    verifDct = dict()
    for annoName, gtfDf in gtfDct.items():

        gtfDf['exonCoord'] = list(zip(gtfDf['start'], gtfDf['end']))

        grpDf = gtfDf.groupby(['seqname', 'strand', 'gene_id', 'transcript_id']).agg({
            'exonCoord': list
        }).reset_index()

        exonVerDf = gtfDf.groupby(
            ['seqname', 'strand', 'gene_id', 'transcript_id']).size().reset_index()
        exonVerDf = exonVerDf.rename(columns={0: 'numExon'})
        exonVerDf = exonVerDf.set_index('transcript_id').drop(
            ['seqname', 'strand', 'gene_id'], axis=1)

        grpDf = grpDf.rename({'exonCoord': 'exonCoordLst'}, axis=1)

        grpDct[annoName] = grpDf
        verifDct[annoName] = exonVerDf

    print()
    print("GTF Grouping Complete!")
    print()

    # Merge all GTFs on transcriptID(and things unique to transcript)
    # Store all possibly different exon coordinates in a new column
    mergeDf = pd.DataFrame()
    mergeSuffix = 1
    for grpDf in grpDct.values():

        if mergeDf.empty:
            mergeDf = grpDf
        else:
            mergeDf = pd.merge(mergeDf, grpDf, how='outer',
                               on=['seqname', 'strand',
                                   'transcript_id', 'gene_id'],
                               indicator='merge_check',
                               suffixes=[mergeSuffix, mergeSuffix+1])

            mergeDf['mergeCheck_' + str(mergeSuffix)] = mergeDf['merge_check']
            mergeDf = mergeDf.drop('merge_check', axis=1)

            mergeSuffix += 1

    print()
    print("GTF Merging Complete!")
    print()

    # Loop through the merge dataframe and compare exon coords across annotations
    rowLst = []
    # 14108
    for row in mergeDf.to_dict('records'):

        multiExonCoordWNanLst = [value for key,
                                 value in row.items() if 'exonCoordLst' in key]

        multiExonCoordLst = [
            exonCoordLst for exonCoordLst in multiExonCoordWNanLst if exonCoordLst is not np.nan]

        if len(multiExonCoordLst) == 1:
            row['outExonLst'] = multiExonCoordLst[0]

        else:

            # Check for a major error: number of exons for the same jxnHash
            # across two annotations is different (this should not happen and
            # would be an issue with id_ujc).
            if any(len(exonCoordLst) != len(multiExonCoordLst[0]) for exonCoordLst in multiExonCoordLst):
                raise Exception("The number of exons for the same jxnHash across "
                                "two annotations is different. jxnHash: "
                                "{}".format(row['transcript_id']))
            else:

                # Otherwise do the following analysis:
                # Loop through each exon and look at the different
                # set of coordinates that exist for each exon
                outExonLst = []
                for exonCoordLst in zip(*multiExonCoordLst):

                    # If this exon has different coords across the species,
                    # find min first start or max last end and set exon to that.
                    if len(set(exonCoordLst)) > 1:

                        startDiff = False
                        endDiff = False

                        if any(exonCoord[0] != exonCoordLst[0][0] for exonCoord in exonCoordLst):

                            startDiff = True

                            minStart = min(exonCoord[0]
                                           for exonCoord in exonCoordLst)
                            endPos = exonCoordLst[0][1]

                            outExonLst.append((minStart, endPos))

                        if any(exonCoord[1] != exonCoordLst[0][1] for exonCoord in exonCoordLst):

                            endDiff = True

                            startPos = exonCoordLst[0][0]
                            maxEnd = max(exonCoord[1]
                                         for exonCoord in exonCoordLst)

                            outExonLst.append((startPos, maxEnd))

                        # Onlt the start of the first exons or the end of the last exons should differ.
                        if startDiff and endDiff:
                            raise Exception(
                                "The start AND end are different for a jxnHash's exon "
                                " across multiple annotation. jxnHash: "
                                "{}".format(row['transcript_id']))

                    else:

                        # If the exon is the same across the different
                        # species -> All is well. Add it to the list.

                        outExonLst.append(exonCoordLst[0])

                row['outExonLst'] = outExonLst

        rowLst.append(row)

    tmpOutDf = pd.DataFrame(rowLst)

    print()
    print("Exon Comparison Across Annotations Complete!")
    print()

    # Used for visual in spyder of how exons are changed
    # kaboom = tmp[tmp['transcript_id'] ==
    #              "000ab5cb90829f133fd4e78f093a2d795bc47c83a7aaef93743560cf7c37f949"].copy()
    # kaboom = kaboom.dropna(axis=1)
    # kaboom = kaboom.explode(
    #     [col for col in kaboom.columns if 'exon' in col or 'Exon' in col])

    # Prepare output GTF
    outDf = tmpOutDf.drop(
        columns=[col for col in tmpOutDf.columns if 'exonCoordLst' in col or 'mergeCheck' in col])

    outDf = outDf.explode('outExonLst')
    outDf[['start', 'end']] = pd.DataFrame(
        outDf['outExonLst'].to_list(), index=outDf.index)

    outDf = outDf[['seqname', 'strand', 'start',
                   'end', 'transcript_id', 'gene_id']]

    # Verify that exons in the input match
    outVerifDf = outDf.groupby(
        ['seqname', 'strand', 'gene_id', 'transcript_id']).size().reset_index()

    outVerifDf = outVerifDf.rename(columns={0: 'numExon'})

    outVerifDf = outVerifDf.set_index('transcript_id').drop(
        ['seqname', 'strand', 'gene_id'], axis=1)

    outVerifDct = outVerifDf.to_dict(orient='index')

    for annoName, exonVerDf in verifDct.items():

        loopDct = exonVerDf.to_dict(orient='index')

        for xscript, numExon in loopDct.items():
            inNumExon = numExon['numExon']
            outNumExon = outVerifDct[xscript]['numExon']

            if inNumExon != outNumExon:
                raise Exception("The number of exons in the output for a "
                                "jxnHash is different. jxnHash: "
                                "{}, annotation: {}".format(xscript, annoName))
    else:
        print("Exon counts are verified!")

    print("Total Number of Unique UJCs in output: {}".format(
        outDf['transcript_id'].nunique()))
    # verify exon counts.

    print("Writing GTF...")
    if os.path.isfile(outGTF):
        os.remove(outGTF)

    trand.io.write_gtf(data=outDf, out_fhs={"gtf": outGTF}, fh_name="gtf")

    print("Script Complete!")

    # Dev code
    # testDf = masterDf[masterDf['transcript_id'] ==
    #                   "000ab5cb90829f133fd4e78f093a2d795bc47c83a7aaef93743560cf7c37f949"]
    # test2Df = gtfDf[gtfDf['transcript_id'] ==
    #                 "000ab5cb90829f133fd4e78f093a2d795bc47c83a7aaef93743560cf7c37f949"]


if __name__ == '__main__':
    global args
    args = getOptions()
    main()
