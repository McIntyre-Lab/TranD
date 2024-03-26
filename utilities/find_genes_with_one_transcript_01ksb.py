#!/usr/bin/env python

import argparse
import trand.io

def getOptions():
    # Parse command line arguments
    parser = argparse.ArgumentParser(description="")
    
    # Input data
    parser.add_argument("-g",
                        "--input-gtf",
                        dest="inGTF", 
                        required=True,
                        help="GTF to find genes with one transcript")
    
    # Output data
    parser.add_argument("-", "--", dest="", required=True, help="")
    
    args = parser.parse_args()
    return args

def main():
    
    inGTF = "/nfshome/k.bankole/mnt/exasmb.rc.ufl.edu-blue/mcintyre/share/sex_specific_splicing/fiveSpecies_annotations/fiveSpecies_2_dmel6_ujc.gtf"
    inDf = trand.io.read_exon_data_from_file(inGTF)
    
    geneDf = inDf.groupby('gene_id').agg({'transcript_id':set}).reset_index()
    
    uniqGeneSet = set(geneDf['gene_id'].to_list())
    
    print("Total number of unique genes: {}".format(len(uniqGeneSet)))
    
    oneXscriptSet = set(geneDf['gene_id'][geneDf['transcript_id'].apply(lambda x: len(x) == 1)].to_list())
    
    print ("Number of genes with one transcript: {}".format(len(set(oneXscriptSet))))    
    print ("Number of genes with more than one transcript (should match number of genes in TranD output): {}".format(len(uniqGeneSet) - len(oneXscriptSet)))
    
    
    # Missing output files...

if __name__ == '__main__':
    # Parse command line arguments
    global args
    args = getOptions()
    main()
