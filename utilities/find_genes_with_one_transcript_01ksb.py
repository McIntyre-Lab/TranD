#!/usr/bin/env python

import argparse
import trand.io
import os
import pandas as pd

def getOptions():
    # Parse command line arguments
    parser = argparse.ArgumentParser(description="Read in GTF. Output list of unique genes (list_{name}_uniq_gene_set) and list of genes with one "
                                     "transcript (list_{name}_genes_w_one_xscript.csv. {name} defualts to name of the input GTF (use --prefix "
                                     "to change it_.")
    
    # Input data
    parser.add_argument("-g",
                        "--input-gtf",
                        dest="inGTF", 
                        required=True,
                        help="GTF to find genes with one transcript")
    
    # Output data
    parser.add_argument("-o",
                        "--output-directory",
                        dest="outDir",
                        required=True,
                        help="Output directory. Must already exist.")
    
    parser.add_argument("-p",
                        "--prefix",
                        dest="prefix",
                        required=False,
                        help="Optional output file prefix")
    
    
    args = parser.parse_args()
    return args

def main():
    
    # inGTF = "/nfshome/k.bankole/mnt/exasmb.rc.ufl.edu-blue/mcintyre/share/sex_specific_splicing/fiveSpecies_annotations/fiveSpecies_2_dmel6_ujc.gtf"
    # outDir = "/nfshome/k.bankole/Desktop"
    
    inGTF = args.inGTF
    outDir = args.outDir
    
    if args.prefix is None:
        prefix = os.path.basename(inGTF).split('.')[0]

    inDf = trand.io.read_exon_data_from_file(inGTF)
    
    
    # Group by gene, count number of unique genes
    geneDf = inDf.groupby('gene_id').agg({'transcript_id':set}).reset_index()
    uniqGeneSet = set(geneDf['gene_id'].to_list())
    print("Total number of unique genes: {}".format(len(uniqGeneSet)))
    
    # Subset for genes with only one transcript
    oneXscriptSet = set(geneDf['gene_id'][geneDf['transcript_id'].apply(lambda x: len(x) == 1)].to_list())
    print ("Number of genes with one transcript: {}".format(len(set(oneXscriptSet))))    
    print ("Number of genes with more than one transcript (should match number of genes in TranD output): {}".format(len(uniqGeneSet) - len(oneXscriptSet)))
    
    # Output both lists
    
    
    outGeneDf = pd.DataFrame(geneDf['gene_id'])
    outGeneDf = outGeneDf.rename({'gene_id':'geneID'})
    
    outOneDf = pd.DataFrame({'geneID':list(oneXscriptSet)}).sort_values(by='geneID').reset_index(drop=True)
    
    geneLstFile = outDir + "/list_{}_uniq_gene_set.csv".format(prefix)
    oneListFile = outDir + "/list_{}_genes_w_one_xscript.csv".format(prefix)
    
    try:
            outGeneDf.to_csv(geneLstFile, index=False)
            outOneDf.to_csv(oneListFile, index=False)
    except OSError:
            raise OSError("Output directory must already exist.")
    
    
    
    
if __name__ == '__main__':
    # Parse command line arguments
    global args
    args = getOptions()
    main()
