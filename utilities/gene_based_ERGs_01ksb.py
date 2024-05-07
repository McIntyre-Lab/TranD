#!/usr/bin/env python

import argparse
import pandas as pd
import trand.io
import re

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

def createExonOutput(forGTF, erInfo):
    
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
    return exonDf


def main():
    
    # eaFile = "/nfshome/k.bankole/mnt/exasmb.rc.ufl.edu-blue/mcintyre/share/sex_specific_splicing/test_gene_based_ERGs/sub_ea_LOC110191367_noIR.csv"
    eaFile = "/nfshome/k.bankole/mnt/exasmb.rc.ufl.edu-blue/mcintyre/share/sex_specific_splicing/test_gene_based_ERGs/sub_keepIR_event_analysis_er.csv"
    # eaFile = "/nfshome/k.bankole/mclab/SHARE/McIntyre_Lab/sex_specific_splicing/trand_1GTF_geneMode_fiveSpecies_ujc/fiveSpecies_2_dser1_ujc_event_analysis_er.csv"
    # pw_eaFile = "/nfshome/k.bankole/mnt/exasmb.rc.ufl.edu-blue/mcintyre/share/sex_specific_splicing/test_gene_based_ERGs/sub_pw_event_analysis.csv"
    # pw_eaFile = "/nfshome/k.bankole/mnt/exasmb.rc.ufl.edu-blue/mcintyre/share/sex_specific_splicing/test_gene_based_ERGs/sub_pw_noIR_event_analysis.csv"
    
    outGTF = "/nfshome/k.bankole/mnt/exasmb.rc.ufl.edu-blue/mcintyre/share/sex_specific_splicing/test_gene_based_ERGs/gene_based_erg_keepIR.gtf"
    # outGTF = "/nfshome/k.bankole/mnt/exasmb.rc.ufl.edu-blue/mcintyre/share/sex_specific_splicing/test_gene_based_ERGs/gene_based_erg_noIR.gtf"
    
    eaFile = args.eaFile
    outGTF = args.outGTF    
    
    
    eaDf = pd.read_csv(eaFile, low_memory=False)       

    
    # pwDf = pd.read_csv(pw_eaFile, low_memory=False)
    # pw_testDf = pwDf[['er_id','er_start']]
    # pw_grp = pw_testDf.groupby('er_id')['er_start'].agg(min)
    
    # Note: can really only ID ERGs for jxnHash stuff, or any transcript IDs that dont contain a | (pipe)
    # Create ERCs for ERGs
    tmpDf = eaDf[['er_id','er_transcript_ids','gene_id']].copy()
    
    testRowLst = [
    {'er_id': '0c10a9919a9cca1fcfd6a40a3e1884251ab221ec924581fdd4ca46a88f951268_exon_1', 'er_transcript_ids': '0c10a9919a9cca1fcfd6a40a3e1884251ab221ec924581fdd4ca46a88f951268', 'gene_id': 'LOC110176799'},
    {'er_id': '0c10a9919a9cca1fcfd6a40a3e1884251ab221ec924581fdd4ca46a88f951268_exon_2', 'er_transcript_ids': '0c10a9919a9cca1fcfd6a40a3e1884251ab221ec924581fdd4ca46a88f951268', 'gene_id': 'LOC110176799'},
    {'er_id': '0c10a9919a9cca1fcfd6a40a3e1884251ab221ec924581fdd4ca46a88f951268_exon_3', 'er_transcript_ids': '0c10a9919a9cca1fcfd6a40a3e1884251ab221ec924581fdd4ca46a88f951268', 'gene_id': 'LOC110176799'},
    {'er_id': '0c10a9919a9cca1fcfd6a40a3e1884251ab221ec924581fdd4ca46a88f951268_exon_4', 'er_transcript_ids': '0c10a9919a9cca1fcfd6a40a3e1884251ab221ec924581fdd4ca46a88f951268', 'gene_id': 'LOC110176799'},
    ]
    
    tmpDf = pd.concat([tmpDf, pd.DataFrame(testRowLst)])
    
    tmpDf.columns = ['ER','transcripts','geneID']
    
    tmpDf['transcripts'] = tmpDf['transcripts'].str.split('|')
    tmpDf = tmpDf.explode('transcripts')
    
    erCatalogue = tmpDf.groupby('transcripts').agg(set).reset_index()
    
    noMultiGeneXscript = erCatalogue['geneID'].apply(lambda x: len(x) == 1).all()
    
    if not noMultiGeneXscript:
        print("There are transcripts belonging to more than one gene. Quitting.")
        quit()
    else:
        erCatalogue['geneID'] = erCatalogue['geneID'].apply(lambda x: list(x)[0])
    
    
    
    erCatalogue['ER'] = erCatalogue.apply(lambda row: {row['geneID'] + ":STE" + x.split('exon_')[1] if 'exon' in x else x for x in row['ER']},
                                          axis =1)

    erCatalogue['ER'] = erCatalogue.apply(lambda row: {x.split(row['geneID'] + ':')[1] for x in row['ER']}, axis=1)
    
    erCatalogue['ER'] = erCatalogue['ER'].apply(lambda x: sorted(x, key=lambda x: int(x.split("ER")[1]) if 'ER' in x else int(x.split("STE")[1])))
    
    erCatalogue['numER'] = erCatalogue['ER'].apply(len)
    
    erCatalogue['ERC'] = erCatalogue.apply(lambda row: 
                                           "{}:{}".format(row['geneID'],
                                                          '_'.join(row['ER'])),
                                           axis=1)
        
    erCatalogue = erCatalogue.sort_values(by=['geneID','numER','ERC'])
        
    ergDf = erCatalogue.groupby('ERC').agg({
        'transcripts':set,
        'geneID':set,
        'numER':set}).reset_index()
    
    noMultiGeneERG = ergDf['geneID'].apply(lambda x: len(x) == 1).all()
    noNumERDiff = ergDf['numER'].apply(lambda x: len(x) == 1).all()
    
    if not noMultiGeneERG:
        print("There are ERGs belonging to more than one gene. Quitting.")
        quit()
    else:
        ergDf['geneID'] = ergDf['geneID'].apply(lambda x: list(x)[0])
    
    if not noNumERDiff:
        print("There are transcripts with differing numbers of ERs in the same group. Quitting.")
        quit()
    else:
        ergDf['numER'] = ergDf['numER'].apply(lambda x: list(x)[0])
    
    ergDf['numTranscripts'] = ergDf['transcripts'].apply(len)
    
    ergDf = ergDf[['ERC','geneID','numER','numTranscripts','transcripts']]
    
    # This probably does not work...    
    row = ergDf.iloc[0]
    
    row['ERC'].split(row['geneID'] + ":")[1]
    int(re.split(r'STE|ER', row['ERC'].split(row['geneID'] + ":")[1]))    
    

        
    a = (2)**(1/3)
    
    ergDf['ERG_id'] = ergDf.apply(lambda row: f'ERG_{row.name + 1}', axis=1)   
    ergDf = ergDf.drop(columns='ERG')
    
    merge = pd.merge(ergDf, erCatalogue, on='ERC', how='outer')
    
    
    
    
    
    # Write GTF
    erInfo = eaDf[['gene_id','er_id','er_chr','er_strand','er_start','er_end']].copy()
    forGTF = merge[['ERG_id', 'ER']].copy()
    forGTF['ER'] = forGTF['ER'].apply(tuple)
    forGTF = forGTF.drop_duplicates()
    trand.io.write_gtf(data=createExonOutput(forGTF=forGTF,erInfo=erInfo), out_fhs={"gtf":outGTF}, fh_name="gtf") 



if __name__ == '__main__':
    # Parse command line arguments
    global args
    args = getOptions()
    main()
