#!/usr/bin/env python

import argparse
import pandas as pd
import trand.io

def getOptions():
        """
        
        Function to store user input via argparse

        Returns
        -------
        args : ARGPARSE ARGUMENTS
                User input via argparse.

        """
        
        # Parse command line arguments
        parser = argparse.ArgumentParser(description="Identifies groups of transcripts that share all of their exon regions.")
        
        # INPUT
        parser.add_argument(
                "-ea",
                "--event-analysis",
                dest="eaFile",
                required=True,
                help="Location of event analysis ER file from gene mode output"
                )
        
        
        # OUTPUT
        parser.add_argument("-o",
                            "--output-gtf",
                            dest="outGTF", 
                            required=True,
                            help="Path and filename for output GTF")
        
        args = parser.parse_args()
        return args

def main():
    
    # eaFile = "/nfshome/k.bankole/mnt/exasmb.rc.ufl.edu-blue/mcintyre/share/sex_specific_splicing/test_gene_based_ERGs/sub_ea_LOC110191367_noIR.csv"
    eaFile = "/nfshome/k.bankole/mnt/exasmb.rc.ufl.edu-blue/mcintyre/share/sex_specific_splicing/test_gene_based_ERGs/sub_keepIR_event_analysis_er.csv"
    # eaFile = "/nfshome/k.bankole/mclab/SHARE/McIntyre_Lab/sex_specific_splicing/trand_1GTF_geneMode_fiveSpecies_ujc/fiveSpecies_2_dser1_ujc_event_analysis_er.csv"
    # pw_eaFile = "/nfshome/k.bankole/mnt/exasmb.rc.ufl.edu-blue/mcintyre/share/sex_specific_splicing/test_gene_based_ERGs/sub_pw_event_analysis.csv"
    # pw_eaFile = "/nfshome/k.bankole/mnt/exasmb.rc.ufl.edu-blue/mcintyre/share/sex_specific_splicing/test_gene_based_ERGs/sub_pw_noIR_event_analysis.csv"
    
    outGTF = "/nfshome/k.bankole/mnt/exasmb.rc.ufl.edu-blue/mcintyre/share/sex_specific_splicing/test_gene_based_ERGs/gene_based_erg_keepIR.gtf"
    outGTF = "/nfshome/k.bankole/mnt/exasmb.rc.ufl.edu-blue/mcintyre/share/sex_specific_splicing/test_gene_based_ERGs/gene_based_erg_noIR.gtf"
    
    eaFile = args.eaFile
    outGTF = args.outGTF
    
    
    
    
    eaDf = pd.read_csv(eaFile, low_memory=False)       

    
    # pwDf = pd.read_csv(pw_eaFile, low_memory=False)
    # pw_testDf = pwDf[['er_id','er_start']]
    # pw_grp = pw_testDf.groupby('er_id')['er_start'].agg(min)
    
    # Note: can really only ID ERGs for jxnHash stuff, or any transcript IDs that dont contain a |


    # Create ERCs for ERGs
    testDf = eaDf[['er_id','er_transcript_ids','gene_id']].copy()
    
    testDf.columns = ['ER','transcripts','geneID']
    
    testDf['transcripts'] = testDf['transcripts'].str.split('|')
    testDf = testDf.explode('transcripts')
    
    erCatalogue = testDf.groupby('transcripts').agg(set).reset_index()
    
    def sortERID(x):
        if "ER" in x:
            return int(x.split(":ER")[1])
        elif "exon" in x:
            return int(x.split("_exon_")[1])
        else:
            return x
        
    erCatalogue['ER'] = erCatalogue['ER'].apply(lambda x: sorted(x,key=lambda x: sortERID(x)))
        
    erCatalogue['ERC'] = erCatalogue['ER'].apply(lambda x: '_'.join(x))
    ergDf = erCatalogue.groupby('ERC')['transcripts'].apply(list).reset_index(name='ERG')
    
    # NEED TO DESPARATELY FIX THIS. ARBITRARY ERG_IDs
    ergDf['ERG_id'] = ergDf.apply(lambda row: f'ERG_{row.name + 1}', axis=1)   
    ergDf = ergDf.drop(columns='ERG')
    
    merge = pd.merge(ergDf, erCatalogue, on='ERC', how='outer')
    
    
    
    
    erInfo = eaDf[['gene_id','er_id','er_chr','er_strand','er_start','er_end']].copy()
    
    forGTF = merge[['ERG_id', 'ER']].copy()
    forGTF['ER'] = forGTF['ER'].apply(tuple)
    forGTF = forGTF.drop_duplicates()
        
    seqnameLst = []
    startLst = []
    endLst = []
    strandLst = []
    ergIDLst = []
    geneIDLst = []
    for row in forGTF.to_dict('records'):
        for er in list(row['ER']):
            seqname = erInfo[erInfo['er_id'] == er]['er_chr'].item()
            strand = erInfo[erInfo['er_id'] == er]['er_strand'].item()
            geneID = erInfo[erInfo['er_id'] == er]['gene_id'].item()
            
            start = erInfo[erInfo['er_id'] == er]['er_start'].item()
            end = erInfo[erInfo['er_id'] == er]['er_end'].item() 
            
            seqnameLst.append(seqname)
            strandLst.append(strand)
            ergIDLst.append(row['ERG_id'])
            geneIDLst.append(geneID)
            
            startLst.append(start)
            endLst.append(end)
                
                        
    exonDf = pd.DataFrame(
            {
                    'seqname':seqnameLst,
                    'start':startLst,
                    'end':endLst,
                    'strand':strandLst,
                    'transcript_id':ergIDLst,
                    'gene_id':geneIDLst
            })
    
    exonDf = exonDf.sort_values(by=['seqname','transcript_id','start']).reset_index(drop=True)

    trand.io.write_gtf(data=exonDf, out_fhs={"gtf":outGTF}, fh_name="gtf")


if __name__ == '__main__':
    # Parse command line arguments
    global args
    args = getOptions()
    main()
