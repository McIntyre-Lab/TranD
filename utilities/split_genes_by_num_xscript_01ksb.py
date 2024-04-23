#!/usr/bin/env python

import argparse
import trand.io
import os
import pandas as pd

def getOptions():
    # Parse command line arguments
    parser = argparse.ArgumentParser(description="Read in GTF. Output list of unique genes below (list_noHeader_lowMem_{name}_gene) and "
                                     "above (list_noHeader_highMem_{name}_gene) an input transcript threshold (--threshold), default of 36. "
                                     "Useful for splitting TranD runs into genes that need low/high dedicated memory."
                                     "{name} defualts to name of the input GTF (use --prefix to change it_.")
    
    # Input data
    parser.add_argument("-g",
                        "--input-gtf",
                        dest="inGTF", 
                        required=True,
                        help="GTF to analyze")
    
    parser.add_argument("-t",
                        "--threshold",
                        dest="threshold", 
                        required=True,
                        default=36,
                        help="Num xscript threshold (<=)")
    
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
    # threshold = int(36)
    # prefix = None
    
    inGTF = args.inGTF
    outDir = args.outDir
    threshold = int(args.threshold)
    prefix = args.prefix
    
    if prefix is None:
        prefix = os.path.basename(inGTF).split('.')[0]

    inDf = trand.io.read_exon_data_from_file(inGTF)
    
    # Group by gene, count number of transcripts per gene
    geneDf = inDf.groupby('gene_id').agg({'transcript_id':'nunique'}).reset_index()
    
    geneDf['transcript_id'] = pd.to_numeric(geneDf['transcript_id'],errors='coerce')
    
    if  geneDf['transcript_id'].isna().any():
        print("something weird happened...")
        quit()
        
        
    uniqGeneSet = set(geneDf['gene_id'].to_list())
    print("Total number of unique genes: {}".format(len(uniqGeneSet)))
    
    lowMemSet = set(geneDf['gene_id'][geneDf['transcript_id'].apply(lambda x: x <= threshold)].to_list())
    highMemSet = set(geneDf['gene_id'][geneDf['transcript_id'].apply(lambda x: x > threshold)].to_list())
    
    if len(lowMemSet) + len(highMemSet) != len(uniqGeneSet):
        print("ERROR: Genes in split lists do not add up to original value. Quitting...")
    else:
        print("Genes successfully split!")
    
    print ("Number of \'low memory\' genes: {}".format(len(set(lowMemSet))))    
    print ("Number of \'high memory\' genes: {}".format(len(set(highMemSet))))    
    
    lowDf = pd.DataFrame({'geneID':list(lowMemSet)}).sort_values(by='geneID').reset_index(drop=True)
    highDf = pd.DataFrame({'geneID':list(highMemSet)}).sort_values(by='geneID').reset_index(drop=True)
    
    lowFile = outDir + "/list_noHeader_lowMem_{}_gene.txt".format(prefix)
    highFile = outDir + "/list_noHeader_highMem_{}_gene.txt".format(prefix)
    
    try:
            lowDf.to_csv(lowFile, header=None, index=False)
            highDf.to_csv(highFile, header=None,index=False)
    except OSError:
            raise OSError("Output directory must already exist.")
    
    
    
    
if __name__ == '__main__':
    # Parse command line arguments
    global args
    args = getOptions()
    main()
