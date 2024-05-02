#!/usr/bin/env python

import argparse
import pandas as pd
import trand.io

def getOptions():
    # Parse command line arguments
    parser = argparse.ArgumentParser(description="Convert event analysis file from TranD 1 GTF gene output "
                                                 "(event_analysis_er.csv) to a GTF where each transcript is "
                                                 "represents all of the exon regions of a gene.")

    # Input data
    parser.add_argument("-e",
                        "--event-analysis",
                        dest="eaFile", 
                        required=True, 
                        help="Path to event analysis file")

    # Output data
    parser.add_argument("-o",
                        "--output-gtf",
                        dest="outGTF", 
                        required=True,
                        help="Path and filename for output GTF")
    
    args = parser.parse_args()
    return args

def main():
    
    # eaFile = "/nfshome/k.bankole/mnt/exasmb.rc.ufl.edu-blue/mcintyre/share/sex_specific_splicing/test_conv_EA2GTF/event_analysis_er.csv"
    # gtfOutPath = "/nfshome/k.bankole/mnt/exasmb.rc.ufl.edu-blue/mcintyre/share/sex_specific_splicing/test_conv_EA2GTF/test_gtf.gtf"
    
    eaFile = args.eaFile
    gtfOutPath = args.outGTF
    
    eaDf = pd.read_csv(eaFile,low_memory=False)
    
    eaDf = eaDf.rename(columns={'er_chr':'seqname',
                        'er_strand':'strand',
                        'er_start':'start',
                        'er_end':'end'})
    
    eaDf['transcript_id'] = eaDf['gene_id']
    eaDf = eaDf[['seqname','start','end','strand','transcript_id','gene_id']]
    outExonDf = eaDf.sort_values(by=['seqname','transcript_id','start']).reset_index(drop=True)

    trand.io.write_gtf(data=outExonDf, out_fhs={"gtf":gtfOutPath}, fh_name="gtf")

if __name__ == '__main__':
    # Parse command line arguments
    global args
    args = getOptions()
    main()
