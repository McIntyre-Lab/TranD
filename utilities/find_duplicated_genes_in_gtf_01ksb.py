#!/usr/bin/env python

import argparse
import pandas as pd
import time
import os
import csv


def getOptions():
    # Parse command line arguments
    parser = argparse.ArgumentParser(description="Analyze a GTF annotation to determine "
                                     "if there are multiple geneIDs with the "
                                     "same start/end. Annotation must have "
                                     "gene features.")

    # Input data
    parser.add_argument("-a",
                        "--annotation",
                        dest="inAnno",
                        required=True,
                        help="Annotation to be analyzed. Must have gene features.")

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

    # drop_columns = ['feature', 'score', 'frame', 'comments']
    # drop_columns = ['source', 'feature', 'score', 'frame', 'comments']
    # drop_columns = ['score', 'frame', 'comments']

    data = pd.read_csv(infile, sep='\t', comment='#',
                       header=None, low_memory=False)
    file_cols = data.columns

    if len(file_cols) < len(all_gtf_columns):
        gtf_cols = all_gtf_columns[:len(file_cols)]
    data.columns = gtf_cols
    # drop_cols = [x for x in drop_columns if x in gtf_cols]

    # data = data[data['feature'] == 'gene']
    # data = data.drop(labels=drop_cols, axis=1)

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

    # missing_value_num = newData.isnull().sum().sum()
    # if missing_value_num > 0:
    #     print("Total number of missing values:", missing_value_num)
    # else:
    #     print("No missing values in data")

    # gene_id_missing_value_num = newData['geneID'].isnull().sum()

    # gene_symbol_missing_value_num = newData['transcriptID'].isnull().sum()

    # if gene_id_missing_value_num > 0:
    #     print("Missing gene_id value number:", gene_id_missing_value_num)
    # if gene_symbol_missing_value_num > 0:
    #     print("Missing transcript_id value number:",
    #           gene_symbol_missing_value_num)

    newData['start'] = pd.to_numeric(newData['start'], downcast="unsigned")
    newData['end'] = pd.to_numeric(newData['end'], downcast="unsigned")

    toc = time.perf_counter()

    print(f"GTF Read complete,  took {toc-omegatic:0.4f} seconds.")
    return newData


def createDupeList(gtfDf):

    # Subset to gene Features
    geneDf = gtfDf[gtfDf['feature'] == "gene"].copy()

    # Keep geneID and attributes together so they are connected during the groupby
    geneDf['geneInfo'] = geneDf.apply(
        lambda x: (x['geneID'], x['attributes']), axis=1)

    # Group by strand and start/end to find genes that have the exact same coordinates
    grpByCoordDf = geneDf.groupby(['strand', 'start', 'end']).agg(
        {'geneInfo': list, 'seqname': set}).reset_index()

    # Check that each set of coordinates is only on one chromosome (just in case some weird GTF is used)
    singleChrDup = grpByCoordDf['seqname'].apply(lambda x: len(x) == 1)
    if not singleChrDup.all():
        raise Exception(
            "There are genes belonging to more than one seqname. Quitting.")
    else:
        grpByCoordDf['seqname'] = grpByCoordDf['seqname'].apply(
            lambda x: list(x)[0])

    # Count number of genes per start/end
    grpByCoordDf['numGenePerCoords'] = grpByCoordDf['geneInfo'].apply(len)

    numGenesWDup = grpByCoordDf[grpByCoordDf['numGenePerCoords']
                                > 1]['numGenePerCoords'].sum()
    pctDup = (numGenesWDup / len(grpByCoordDf))

    print(f"{numGenesWDup} ({pctDup:.4%} of total) have duplicate geneIDs (2+ genes with the same coordinates and multiple geneIDs)")

    # Subset to only genes with dupes
    dupeDf = grpByCoordDf[grpByCoordDf['numGenePerCoords']
                          > 1].drop('numGenePerCoords', axis=1)

    # Raise a custom error to stop the program if there are no duplicate genes
    if len(dupeDf) == 0:
        raise ValueError("NO DUPE ERROR")

    # Sort geneIDs alphabetically
    dupeDf['geneInfo'] = dupeDf['geneInfo'].apply(
        lambda x: sorted(x, key=lambda y: y[0]))

    dupeDf = dupeDf.rename(columns={'geneInfo': 'geneInfo_ORIG'})

    # For dupes, pick the first geneID (and related attributes) alphabetically
    # to represent all genes with the exact same coordinates
    dupeDf['geneID_NEW'] = dupeDf['geneInfo_ORIG'].apply(lambda x: x[0][0])
    dupeDf['attributes_NEW'] = dupeDf['geneInfo_ORIG'].apply(lambda x: x[0][1])

    # Make df unique on original geneID, split "info" column into geneID and attributes
    dupeDf = dupeDf.explode('geneInfo_ORIG')
    dupeDf[['geneID_ORIG', 'attributes_ORIG']] = pd.DataFrame(
        dupeDf['geneInfo_ORIG'].tolist(), index=dupeDf.index)

    # Reorganize and cleanup columns

    dupeDf = dupeDf[['geneID_ORIG', 'attributes_ORIG', 'geneID_NEW',
                     'attributes_NEW', 'start', 'end', 'strand', 'seqname']]

    return dupeDf


def main():

    inAnno = "/TB14/TB14/blue_copy/references/dmel_fb650/dmel-all-r6.50.gtf"
    # inAnno = "/TB14/TB14/blue_copy/references/dsim_fb202/dsim-all-r2.02.gtf"
    # inAnno = "/TB14/TB14/blue_copy/references/dser1.1/GCF_002093755.2/genomic.gtf"
    # inAnno = "/TB14/TB14/blue_copy/references/dsan_Prin_1.1/GCF_016746245.2/genomic.gtf"
    # inAnno = "/TB14/TB14/blue_copy/references/dyak_Prin_Tai18E2_2.1/GCF_016746365.2/genomic.gtf"
    prefix = None
    outdir = "/nfshome/k.bankole/Desktop/test_folder"

    inAnno = args.inAnno
    prefix = args.prefix
    outdir = args.outdir

    # Read GTF
    gtfDf = read_all_gtf_data_from_file(inAnno)

    print("Num unique genes:", gtfDf['geneID'].nunique())

    # Create list of duplicate genes. If there are no duplicate genes print and quit.
    try:
        dupeDf = createDupeList(gtfDf=gtfDf)
    except ValueError as error:
        if str(error) == "NO DUPE ERROR":
            print(
                "\nThere are no duplicate genes in this GTF! No output will be created.")
        else:
            raise

    # Create a list of genes and their xscripts
    gene2XscriptDf = gtfDf[['transcriptID', 'geneID', 'attributes']
                           ].drop_duplicates().dropna(ignore_index=True).copy()

    # Subset this to the genes that have dupes -> store a list of transcripts
    # that are from duplicate genes (unique on transcriptID)
    dupeXscriptDf = gene2XscriptDf[gene2XscriptDf['geneID'].isin(
        dupeDf['geneID_ORIG'])]

    # 'Fix' original GTF by replacing duplicate genes with the first alphabetical selected earlier
    rowLst = []
    for row in gtfDf.to_dict('records'):

        if row['geneID'] in dupeDf['geneID_ORIG'].tolist():

            replaceInfoDct = dupeDf[dupeDf['geneID_ORIG']
                                    == row['geneID']].squeeze().to_dict()

            # Check that seqname and strand are the same?
            # I'm not sure why they wouldnt be but it doesn't hurt to check

            if replaceInfoDct['seqname'] != row['seqname'] and replaceInfoDct['strand'] != row['strand']:
                raise Exception("For some reason, though they have the same coordinates, "
                                f"the seqname/strand for {row['geneID']} in the original GTF "
                                f"does not match the seqname/strand for {replaceInfoDct['geneID_NEW']}. Quitting.")

            row['geneID'] = replaceInfoDct['geneID_NEW']
            row['attributes'] = replaceInfoDct['attributes_NEW']

        rowLst.append(row)

    # New GTF using 'fixed' geneIDs
    newGTFDf = pd.DataFrame(rowLst)

    # Check that this new GTF has the correct number of geneIDs:
    # number of unique genes originally - (number of genes with duplicates - number of genes used to represent the duplicates)

    removedGenes = set(dupeDf['geneID_ORIG'].tolist()) - \
        set(dupeDf['geneID_NEW'].tolist())

    correctNewNumGene = gtfDf['geneID'].nunique() - len(removedGenes)
    if newGTFDf['geneID'].nunique() != correctNewNumGene:
        raise Exception(
            "Something went wrong. The output GTF does not have the correct number of genes.")
    else:

        # Drop duplicates to remove duplicate gene features
        outGTFDf = newGTFDf.drop_duplicates().copy()

    # If the removed genes are still in any of the GTF rows throw an error
    if any(outGTFDf.isin(removedGenes).any()):
        raise Exception("Dupes were not correctly removed. Quitting.")

    # Create attributes column for output
    # Create a quick function to put transcriptID and geneID back in attributes
    # (for readability I promise)
    def createAttributes(x):
        if x['transcriptID']:
            return "transcript_id \"{}\"; gene_id \"{}\"; {}".format(x['transcriptID'], x['geneID'], x['attributes'])
        else:
            return "gene_id \"{}\"; {}".format(x['geneID'], x['attributes'])
    outGTFDf['attribute'] = outGTFDf.apply(
        lambda x: createAttributes(x), axis=1)

    # OUTPUT PREFIX
    if prefix:
        outPrefix = "{}/{}".format(outdir, prefix)
    else:
        fileName = os.path.basename(inAnno).split('.gtf')[0]
        outPrefix = "{}/{}".format(outdir, fileName)

    # OUTPUT CSVs
    dupeOutFile = f"{outPrefix}_duplicated_genes.csv"
    dupeDf.to_csv(dupeOutFile, index=False)

    dupXscrOutFile = f"{outPrefix}_transcripts_of_duplicated_genes.csv"
    dupeXscriptDf.to_csv(dupXscrOutFile, index=False)

    # OUTPUT NEW GTF
    outColLst = ['seqname', 'source', 'feature', 'start', 'end', 'score', 'strand',
                 'frame', 'attribute']
    outGTFDf = outGTFDf.reindex(columns=outColLst)

    gtfOutFile = f"{outPrefix}_noDup.gtf"

    outGTFDf.to_csv(gtfOutFile, sep="\t", index=False, header=False,
                    doublequote=False, quoting=csv.QUOTE_NONE)


if __name__ == '__main__':
    # Parse command line arguments
    global args
    args = getOptions()
    main()
