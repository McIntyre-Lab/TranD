#!/usr/bin/env python

import argparse
import pandas as pd
import time


def getOptions():
    # Parse command line arguments
    parser = argparse.ArgumentParser(description="Outputs a list of genes that have a mismatch "
                                     "event (a transcript switches genes from "
                                     "original after GFFCompare) post-GFFCompare "
                                     "for a self-mapped GTF (ex: dmel650_2_dmel6). "
                                     "Also prints how many transcripts/genes have "
                                     "a mismatch event")

    # Input data
    parser.add_argument("-a",
                        "--annotation",
                        dest="annoFile",
                        required=True,
                        help="Reference annotation input to GFFCompare.")

    parser.add_argument("-pc",
                        "--protein-coding",
                        dest="pcGTFFile",
                        required=True,
                        help="Refrence annotation subset to only protein coding")

    parser.add_argument("-k",
                        "--gene-key",
                        dest="geneKeyFile",
                        required=True,
                        help="Gene Key output from GFFCompare")
    # Output data
    parser.add_argument("-p",
                        "--prefix",
                        dest="prefix",
                        required=False,
                        help="Name for output files. Optional. Defaults to name of original file.")

    parser.add_argument("-o",
                        "--outdir",
                        dest="outdir",
                        required=True,
                        help="Output directory. Must already exist.")

    args = parser.parse_args()
    return args


def read_all_gtf_data_from_file(infile):
    """
    Create a pandas dataframe of a GTF file creating columns specific for geneID and transcriptID

    Based on trand.io.read_exon_data_from_file
    """

    print("Reading GTF...")

    omegatic = time.perf_counter()

    all_gtf_columns = ['seqname', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame',
                       'attributes', 'comments']

    data = pd.read_csv(infile, sep='\t', comment='#',
                       header=None, low_memory=False)
    file_cols = data.columns

    if len(file_cols) < len(all_gtf_columns):
        gtf_cols = all_gtf_columns[:len(file_cols)]
    data.columns = gtf_cols

    data['seqname'] = data['seqname'].astype(str)
    data['source'] = data['source'].astype(str)
    data['feature'] = data['feature'].astype(str)
    data['start'] = data['start'].astype(int)
    data['end'] = data['end'].astype(int)
    data['attributes'] = data['attributes'].astype(str)

    data.reset_index(drop=True, inplace=True)

    seqnameLst = []
    sourceLst = []
    featureLst = []
    startLst = []
    endLst = []
    scoreLst = []
    strandLst = []
    frameLst = []
    geneIDLst = []
    transcriptIDLst = []
    otherAttrLst = []

    for row in data.to_dict('records'):

        rawAttr = row['attributes']
        attrLst = [x.strip() for x in rawAttr.strip().split(';')]

        gnTrAttr = [
            x for x in attrLst if 'gene_id' in x or 'transcript_id' in x]
        otherAttr = [x for x in attrLst if x and x not in gnTrAttr]

        gene_id, transcript_id = None, None

        for item in gnTrAttr:

            if 'gene_id' in item:
                gene_id = item.split('gene_id')[1].strip().strip('\"')

            elif 'transcript_id' in item:
                transcript_id = item.split(
                    'transcript_id')[-1].strip().strip('\"')

        if not gene_id:
            print("gene_id not found in:", row)
            gene_id = None

        if not transcript_id and row['feature'] != 'gene':
            print("gene_symbol not found in:", row)
            transcript_id = None

        seqnameLst.append(row['seqname'])
        sourceLst.append(row['source'])
        featureLst.append(row['feature'])
        startLst.append(row['start'])
        endLst.append(row['end'])
        scoreLst.append(row['score'])
        strandLst.append(row['strand'])
        frameLst.append(row['frame'])

        geneIDLst.append(gene_id)
        transcriptIDLst.append(transcript_id)
        otherAttrLst.append("; ".join(otherAttr))

    newData = pd.DataFrame(
        {
            'seqname': seqnameLst,
            'source': sourceLst,
            'feature': featureLst,
            'start': startLst,
            'end': endLst,
            'score': scoreLst,
            'strand': strandLst,
            'frame': frameLst,
            'geneID': geneIDLst,
            'transcriptID': transcriptIDLst,
            'attributes': otherAttrLst
        })

    print("GTF rows:", newData.shape[0])

    newData['start'] = pd.to_numeric(newData['start'], downcast="unsigned")
    newData['end'] = pd.to_numeric(newData['end'], downcast="unsigned")

    toc = time.perf_counter()

    print(f"GTF Read complete,  took {toc-omegatic:0.4f} seconds.")
    return newData


def main():

    annoFile = "/TB14/TB14/blue_copy/references/dmel_fb650/dmel650_noDup.gtf"
    pcGTFFile = "/TB14/TB14/blue_copy/references/dmel_fb650/dmel650_protein_coding_ref.gtf"
    geneKeyFile = "/TB14/TB14/blue_copy/references/dmel_fb650/dmel650_2_dmel6_testNoDup_gffcompare_gene_key.csv"

    # annoFile = args.annoFile
    # pcGTFFile = args.pcGTFFile
    # geneKeyFile = args.geneKeyFile

    # Read in annotation and protein coding annotation
    annoDf = read_all_gtf_data_from_file(annoFile)
    pcGTFDf = read_all_gtf_data_from_file(pcGTFFile)

    # Make a DF that is just every gene and transcript pairing
    xscript2GeneDf = annoDf[['geneID', 'transcriptID']].dropna()
    xscript2GeneDf = xscript2GeneDf.rename(
        columns={'geneID': 'geneID_ORIG'}).drop_duplicates().reset_index(drop=True)

    # Make a list of unique protein coding genes using protein coding GTF
    proteinCodingGeneLst = pcGTFDf['geneID'].drop_duplicates().tolist()

    # flag protein coding genes using above list
    xscript2GeneDf['flag_proteinCoding'] = xscript2GeneDf['geneID_ORIG'].isin(
        proteinCodingGeneLst).astype(int)

    # Read in GFFCompare Gene Key and change it to transcriptID and *output* geneID
    keyDf = pd.read_csv(geneKeyFile, low_memory=False)
    keyDf = keyDf[['transcript_id', 'output_gene_id']]
    keyDf = keyDf.rename(
        columns={'transcript_id': 'transcriptID', 'output_gene_id': 'geneID_postGffc'})

    # Merge on transcript ID to compare genes pre- and post- gffcompare
    mergeDf = pd.merge(keyDf, xscript2GeneDf, on='transcriptID',
                       how='outer', indicator='merge_check')

    mergeDf = mergeDf[['transcriptID', 'geneID_ORIG',
                       'geneID_postGffc', 'flag_proteinCoding', 'merge_check']]

    # Flag mismatches
    mergeDf['flag_geneMismatch'] = (
        mergeDf['geneID_postGffc'] != mergeDf['geneID_ORIG']).astype(int)

    # there should be no left_only. right_only=unmapped (this is verified)
    if any(mergeDf['merge_check'] == 'left_only'):
        raise Exception("An error occurred during the merge.")

    # Subset to both
    bothDf = mergeDf[mergeDf['merge_check']
                     == 'both'].drop('merge_check', axis=1)

    # Count num transcripts diff post-GFFC
    numTrNonMatch = bothDf[bothDf['flag_geneMismatch'] == 1].count()[
        'flag_geneMismatch']
    numTrMatch = bothDf[bothDf['flag_geneMismatch'] == 0].count()[
        'flag_geneMismatch']

    pctTrMismatch = numTrNonMatch / (numTrMatch + numTrNonMatch)

    print(f"There are {numTrNonMatch} ({pctTrMismatch:0.4%} of total) transcripts that "
          "have different genes after GFFCompare.")

    # Count num genes diff post-GFFC

    bothDf['geneID_ORIG'].nunique()
    17449

    geneOnlyDf = bothDf.groupby('')

    geneOnlyDf = bothDf.drop('transcriptID', axis=1).drop_duplicates()

    numGnNonMatch = geneOnlyDf[geneOnlyDf['flag_geneMismatch']
                               == 1]['geneID_ORIG'].nunique()

    numGnMatch = geneOnlyDf[geneOnlyDf['flag_geneMismatch']
                            == 0]['geneID_ORIG'].nunique()

    pctGnMismatch = numGnNonMatch / (numGnMatch + numGnNonMatch)

    print(f"There are {numGnNonMatch} ({pctGnMismatch:0.4%} of total) genes involved "
          "in these mismatch events. ")


if __name__ == '__main__':
    # Parse command line arguments
    global args
    args = getOptions()
    main()


# Below is a list of genes that have the same coordinates as another gene on the list, but on different strands
# Looking at these shows that GFFCompare will correctly track strand; no transcript is flipped to a gene on another strand
# test = mergeDf[mergeDf['geneID_NODUP_ORIG'].isin([
#     'FBgn0031758',
#     'FBgn0263442',
#     'FBgn0262819',
#     'FBgn0265363',
#     'FBgn0262290',
#     'FBgn0262462',
#     'FBgn0263409',
#     'FBgn0267612'
# ])]

# left_only = not an exon. right_only = unmapped
