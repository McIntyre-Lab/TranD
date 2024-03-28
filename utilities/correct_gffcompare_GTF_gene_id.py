#!/usr/bin/env python

import argparse
import pandas as pd
import numpy as np
import csv
import sys

"""

3/27: KSB updated so that it does not quit upon trans-spliced genes, just removes them from the xcrptDF (GFFCompare annotated GTF) and bins them in separate output.
        - requires -b paramater to do above

"""

def getOptions():
    # Parse command line arguments
    parser = argparse.ArgumentParser(description="Correct the gene_id values of the GTF-like output from GFFcompare (v12.2 or higher) to match the gene_id values associated with the transcripts")

    # Input data
    parser.add_argument(
        "-a",
        "--annotation",
        dest="inAnnot",
        required=True,
        help=(
            "Input reference annotated GTF-like file from GFFcompare >= v12.2 (*.annotated.gtf). "
            "NOTE: If multiple GTF files need to have gene_id values assigned, must provide "
            "gffcompare with a single concatenated GTF where each individual file has a sampleID "
            "in the \'source\' column (column 2) of each GTF file."
        )
    )
    parser.add_argument(
        "-k",
        "--keep-names",
        dest="inKeep",
        required=False,
        action="store_true",
        help=(
            "When a reference gene_id was not associated, use original gene_id "
            "(Default: Replace with GFFcompare XLOC values)."
        )
    )
    parser.add_argument(
        "-n",
        "--use-gene_name",
        dest="inName",
        required=False,
        action="store_true",
        help=(
            "When a ref_gene_id is not present in GFFcompare output for class "
            "codes not in [s, x, y, i, p, r, u], check for ref_gene_name or "
            "gene_name to use instead (Default: Only use ref_gene_id, else "
            "replace with GFFcompare XLOC values or original gene_id if -k is used)."
        )
    )
    parser.add_argument(
        "-g",
        "--gtf",
        dest="inGTF",
        required=True,
        help=(
            "Input GTF file to be modified (original GTF given to GFFcompare). "
            "NOTE: The name of this file MUST match what was given to "
            "gffcompare so the correct transcripts can be identified in the "
            "output if multiple files were provided to gffcompare."
        )
    )

    # Output data
    parser.add_argument(
        "-o",
        "--output",
        dest="outFile",
        required=True,
        help="Output GTF file name for isoform exons with associated gene_id (file is UNSORTED)"
    )
    parser.add_argument(
        "--key",
        dest="outKey",
        required=True,
        help=(
            "Output CSV key file for each transcript_id, input_gene_id, and "
            "output_gene_id, where transcript_id is the list of all unique "
            "transcripts in the output GTF file, input_gene_id is the gene_id "
            "of the transcript in the input GTF file, and output_gene_id is the "
            "new associated gene_id based on the gffcompare annotation file and "
            "associated with the transcript_id in the output GTF file."
        )
    )
    
    parser.add_argument(
        "--bin-TS",
        dest="binTS",
        required=False,
        default=False,
        help=("Add a path for binned trans-spliced transcripts. will store trans-spliced transcripts separately, rather than stopping"
              " the whole script (default behavior)."))

    args = parser.parse_args()
    return args

def main():

    # inAnnot = "~/mnt/exasmb.rc.ufl.edu-blue/mcintyre/share/transcript_ortholog/test/test.annotated.gtf"
    # # inAnnot = "~/mnt/exasmb.rc.ufl.edu-blue/mcintyre/share/transcript_ortholog/test/test.gtf"
    # keepOrigGene = False
    # useGeneName = False
    # inGTF = "/nfshome/k.bankole/mnt/exasmb.rc.ufl.edu-blue/mcintyre/share/transcript_ortholog/yak_2_dyak2_uniq_jxnHash_noGeneID.gtf"
    # binTS = True
    
    inAnnot = args.inAnnot
    keepOrigGene = args.inKeep
    useGeneName = args.inName
    inGTF = args.inGTF
    binTS = args.binTS
    
    # Get gffcompare annotation output 
    annotDF = pd.read_csv(inAnnot,names=['chr','source','feature','start','end','score',
                                        'strand','frame','attribute'], sep="\t",low_memory=False)
    print("Total lines in gffcompare output = "+str(len(annotDF)))
    print("Total transcripts in gffcompare output = "+str(len(annotDF[annotDF["feature"]=="transcript"])))
    print("Total input sources contributing to gffcompare output = "+str(annotDF["source"].nunique()))

    # Get gene_id and transcript_id values from attributes column
    # Get transcript_id to gene_id pairs using transcript features of gffcompare output
    # If ref_gene_id attribute present (from gffcompare v12.2), assign as gene_id
    # else if xloc attribute present, assign as gene_id (if requested)
    # else not properly formatted
    # Updated to be faster. Verified it has the same behavior as previously.
    
    xcrptDF = annotDF[annotDF["feature"]=="transcript"].reset_index(drop=True)
    
    chrLst = []
    sourceLst = []
    featureLst = []
    startLst = []
    endLst = []
    scoreLst = []
    strandLst = []
    frameLst = []
    attributeLst = []
    geneIDLst = []
    xscriptIDLst = []
    
    for row in xcrptDF.to_dict('records'):
        rawAttr = row['attribute']
        attrLst = [x.strip() for x in rawAttr.strip().split(';')] 
        ## previous gffcompare versions do not have "ref_gene"
        exp_attrs = [x for x in attrLst if 'transcript_id' in x or 'gene_id' in x or 'xloc' in x or "gene_name" in x or "class_code" in x]
        transcript_id, orig_gene_id, ref_gene_id, ref_gene_name, gene_name, xloc, classCode = np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan
        for item in exp_attrs:
            if 'transcript_id' in item:
                transcript_id = item.split('transcript_id')[-1].strip().strip('\"')
            elif 'ref_gene_id' in item:
                ref_gene_id = item.split('ref_gene_id')[-1].strip().strip('\"')
            elif 'ref_gene_name' in item:
                ref_gene_name = item.split('ref_gene_name')[-1].strip().strip('\"')
            elif 'gene_name' in item:
                gene_name = item.split('gene_name')[-1].strip().strip('\"')
            elif 'gene_id' in item:
                orig_gene_id = item.split('gene_id')[-1].strip().strip('\"')
            elif 'xloc' in item:
                xloc = item.split('xloc')[-1].strip().strip('\"')
            elif 'class_code' in item:
                classCode = item.split('class_code')[-1].strip().strip('\"')
        if transcript_id is np.nan:
            print("WARNING: transcript_id not found in {}".format(row))
        
        if ref_gene_id is np.nan:
            if not useGeneName:
                if not keepOrigGene:
                    if xloc is np.nan:
                        print("WARNING: no attribute for gene_id assignment (ref_gene_id or xloc) found in {}".format(row))
                    else:
                        gene_id = xloc
                else:
                    gene_id = orig_gene_id
            else:
                if classCode not in ['s', 'x', 'y', 'i', 'p', 'r', 'u']:
                    if ref_gene_name is np.nan and gene_name is np.nan:
                        if not keepOrigGene:
                            if xloc is np.nan:
                                print("WARNING: no attribute for gene_id assignment (ref_gene_id or xloc) found in {}".format(row))
                            else:
                                gene_id = xloc
                        else:
                            gene_id = orig_gene_id
                    else:
                        if ref_gene_name is np.nan:
                            gene_id = gene_name
                        else:
                            gene_id = ref_gene_name
        else:
            gene_id = ref_gene_id
            
        
        chrLst.append(row['chr'])
        sourceLst.append(row['source'])
        featureLst.append(row['feature'])
        startLst.append(row['start'])
        endLst.append(row['end'])
        scoreLst.append(row['score'])
        strandLst.append(row['strand'])
        frameLst.append(row['frame'])
        attributeLst.append(row['attribute'])
        
        geneIDLst.append(gene_id)
        xscriptIDLst.append(transcript_id)
        
    xcrptDF = pd.DataFrame({
        'chr': chrLst,
        'source': sourceLst,
        'feature': featureLst,
        'start': startLst,
        'end': endLst,
        'score': scoreLst,
        'strand': strandLst,
        'frame': frameLst,
        'attribute': attributeLst,
        'gene_id': geneIDLst,
        'transcript_id': xscriptIDLst
    })

    # Get GTF file and count total lines
    gtfDf = pd.read_csv(inGTF,names=['chr','source','feature','start','end','score',
                                        'strand','frame','attribute'], sep="\t",low_memory=False,comment="#")
    print("Total lines in original GTF= "+str(len(gtfDf)))
    
    
    # Get gene_id and transcript_id values from attributes column
    # An updated method that is much faster . Verified it has the same behavior as the old method.
    
    seqnameLst = []
    sourceLst = []
    featureLst = []
    startLst = []
    endLst = []
    scoreLst = []
    strandLst = []
    frameLst = []
    attributeLst = []
    geneIDLst = []
    gIdxLst = []
    xscriptIDLst = []
    tIdxLst = []
    
    for row in gtfDf.to_dict('records'):
            rawAttr = row['attribute']
            attrLst = [x.strip() for x in rawAttr.strip().split(';')]
            gnTrAttr = [x for x in attrLst if 'transcript_id' in x or 'gene_id' in x]
            gene_id, transcript_id = None, None
            
            for item in gnTrAttr:
                    if 'gene_id' in item:
                            gene_id = item.split('gene_id')[1].strip().strip('\"')
                            gene_idx = gnTrAttr.index(item)

                    elif 'transcript_id' in item:
                            transcript_id = item.split('transcript_id')[-1].strip().strip('\"')
                            transcript_idx = gnTrAttr.index(item)

            if not gene_id:
                    print("gene_id not found in '{}'", row)
                    gene_id = None
                    
            if not transcript_id:
                    print("transcript_id not found in '{}'", row)
                    transcript_id = None

            
            seqnameLst.append(row['chr'])
            sourceLst.append(row['source'])
            featureLst.append(row['feature'])
            startLst.append(row['start'])
            endLst.append(row['end'])
            scoreLst.append(row['score'])
            strandLst.append(row['strand'])
            frameLst.append(row['frame'])
            attributeLst.append(row['attribute'])
            
            geneIDLst.append(gene_id)
            gIdxLst.append(gene_idx)
            xscriptIDLst.append(transcript_id)
            tIdxLst.append(transcript_idx)

    
    newGTFDf = pd.DataFrame(
            {
                    'chr':seqnameLst,
                    'source':sourceLst,
                    'feature':featureLst,
                    'start':startLst,
                    'end':endLst,
                    'score':scoreLst,
                    'strand':strandLst,
                    'frame':frameLst,
                    'attribute':attributeLst,
                    'gene_id':geneIDLst,
                    'gene_idx':gIdxLst,
                    'transcript_id':xscriptIDLst,
                    'transcript_idx':tIdxLst
            })
            
    print("Total transcripts in original GTF= "+str(newGTFDf["transcript_id"].nunique()))
    print("Total genes in original GTF= "+str(newGTFDf["gene_id"].nunique()))
    

    # Check for unique pairs of gene and transcript
    # Duplicate values could be due to trans-spliced inputs
    if xcrptDF["transcript_id"].duplicated().any():
        
        print("!!!ERROR: Duplicate gene_id:transcript_id pairs in annotation file"
              "...can be due to trans-spliced inputs.\n"
              "\tThe following are the rows with duplicated gene-transcript:\n{}".format(
                  xcrptDF[xcrptDF["transcript_id"].duplicated(keep=False)].to_string()))
        
        if binTS:
            dropDupeDf = xcrptDF.drop_duplicates(subset=['transcript_id'], keep=False)
            droppedRowDf = xcrptDF[~xcrptDF.index.isin(dropDupeDf.index)]
            
            print("Number of removed trans-spliced transcripts: ", droppedRowDf['transcript_id'].nunique())
            
            xcrptDF = dropDupeDf.copy(deep=True).reset_index(drop=True)
            
            
            # Remove dropped transcripts from the GTF as well
            dropXscriptLst = droppedRowDf['transcript_id'].unique().tolist()
            newGTFDf = newGTFDf[~newGTFDf['transcript_id'].isin(dropXscriptLst)]
            
            droppedRowDf.to_csv(binTS,index=False)
            
        else:
            exit()
    
    # Merge gene_id from gffcompare output to original gtf on transcript_id
    mergeDFxcrpt = pd.merge(
            newGTFDf,
            xcrptDF[["transcript_id", "gene_id"]],
            how="left",
            on="transcript_id",
            suffixes=["_gtf", "_cmp"],
            indicator="merge_check",
            validate="m:1"
    )
    if mergeDFxcrpt["merge_check"].value_counts()["left_only"] > 0:
        print("WARNING: Original GTF contains transcripts not present in gffcompare output.")
    
    # Check for proper merge
    if len(mergeDFxcrpt[mergeDFxcrpt['gene_id_cmp'].isna()]) > 0:
        sys.exit("!!! ERROR : UNEXPECTED MERGE RESULTING IN MISSING gene_id ASSOCIATION")
    
    
    
    
    # New updated faster version. Verified it has the same behavior as before
    # Create new attribute column with associated_gene as gene_id

    chrLst = []
    sourceLst = []
    featureLst = []
    startLst = []
    endLst = []
    scoreLst = []
    strandLst = []
    frameLst = []
    attributeLst = []
    
    for row in mergeDFxcrpt.to_dict('records'):
        
        if row['gene_idx'] == 0:
            new_attribute = "gene_id \"" + row['gene_id_cmp'] + "\";" + ";".join(row['attribute'].split(";")[1:])
        else:
            new_attribute = ";".join(row['attribute'].split(";")[0:int(row['gene_idx'])]) + "; gene_id \"" + row['gene_id_cmp'] + "\";" + ";".join(row['attribute'].split(";")[int(row['gene_idx'])+1:])
              
        chrLst.append(row['chr'])
        sourceLst.append(row['source'])
        featureLst.append(row['feature'])
        startLst.append(row['start'])
        endLst.append(row['end'])
        scoreLst.append(row['score'])
        strandLst.append(row['strand'])
        frameLst.append(row['frame'])
        attributeLst.append(new_attribute)
    
    outGTFDf = pd.DataFrame({
        'chr': chrLst,
        'source': sourceLst,
        'feature': featureLst,
        'start': startLst,
        'end': endLst,
        'score': scoreLst,
        'strand': strandLst,
        'frame': frameLst,
        'new_attribute': attributeLst
    })
    
    
    
    outGTFDf[['chr','source','feature','start','end','score','strand','frame',
             'new_attribute']].to_csv(args.outFile,sep="\t",index=False,header=False,
             doublequote=False,quoting=csv.QUOTE_NONE)

    # Output key file for transcript_id to input_gene_id to output_gene_id
    keyDF = mergeDFxcrpt[["transcript_id", "gene_id_gtf", "gene_id_cmp"]].drop_duplicates().rename(
        columns={"gene_id_gtf": "input_gene_id","gene_id_cmp": "output_gene_id"})
    keyDF.to_csv(args.outKey, index=False)


if __name__ == '__main__':
    # Parse command line arguments
    global args
    args = getOptions()
    main()

