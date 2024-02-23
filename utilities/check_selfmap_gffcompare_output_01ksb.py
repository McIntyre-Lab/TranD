#!/usr/bin/env python

import argparse
import trand.io
import pandas as pd
import numpy as np

def getOptions():
        # Parse command line arguments
        parser = argparse.ArgumentParser(description="Creates a summary file for the differences between the geneIDs "
                                         "pre- and post- gffcompare. Takes in: 1. "
                                         "Takes in: 1. The gene key from correct gffcompare output "
                                         "2. Reference annotation (used when doing gffcompare) "
                                         "3. A GTF with all transcripts that did not map. "
                                         "Using an existing output directory and file prefix, "
                                         "Outputs: 1. A file listing all the transcripts that should be removed from the corrected_associated GTF - "
                                         "prefix_remove_these_transcripts_from_corrected_gtf.txt "
                                         "2. A file that lists a reference transcript, its reference geneID, its geneID in the corrected output "
                                         ", and a whether it succsesfully mapped. Also includes a flag that shows if gffcompare changed "
                                         "the transcripts geneID. - prefix_gffcompare_correction_summary.csv"
                                         "Situations where a transcript is to be removed from the corrected GTF: "
                                         "1. The transcripts are not unmapped, but only appear in the reference annotation.")        
        # Input data
        parser.add_argument("-k",
                            "--gene-key",
                            dest="inGeneKey",
                            required=True, 
                            help="Full path to gene key from correct_gffcompare output")
        
        
        parser.add_argument("-r",
                            "--reference-annotation",
                            dest="inRefAnno",
                            required=True, 
                            help="Full path to reference GTF (that was used as the reference for GFFCompare)")
        
        parser.add_argument("-un",
                            "--unmapped-gtf",
                            dest="inUnmapGTF",
                            required=True, 
                            help="Full path GTF containing information on transcripts that did not successfully map to genome")
        
        # Output data
        parser.add_argument("-p",
                            "--filePrefix", 
                            dest="filePrefix", 
                            required=True,
                            help="File prefix for output files.")
        
        parser.add_argument("-o",
                            "--outDir", 
                            dest="outDir", 
                            required=True,
                            help="Path to output directory (must already exist).")
        
        
        args = parser.parse_args()
        return args

def main():
        inGeneKey = args.inGeneKey
        inRefAnno = args.inRefAnno
        inUnmapGTF = args.inUnmapGTF
        outDir = args.outDir
        filePrefix = args.filePrefix
        
        # species = "ser"
        # species = "sim"
        # species = "yak"
        # species = "san"
        # species = "mel"

        # if species == "mel":        
        #         inGeneKey = "/nfshome/k.bankole/mnt/exasmb.rc.ufl.edu-blue/mcintyre/share/references/dmel_fb650/dmel650_2_dmel6_gene_key.csv"
        #         inRefAnno = "/nfshome/k.bankole/mnt/exasmb.rc.ufl.edu-blue/mcintyre/share/references/dmel_fb650/dmel-all-r6.50.gtf"
        #         inUJCIndex = "/nfshome/k.bankole/mnt/exasmb.rc.ufl.edu-blue/mcintyre/share/references/dmel_fb650/dmel650_2_dmel6_ujc_xscript_index.csv"
        #         inUnmapGTF = "/nfshome/k.bankole/mnt/exasmb.rc.ufl.edu-blue/mcintyre/share/references/dmel_fb650/dmel650_2_dmel6_unmapped.gtf"
        
        # elif species == "san":
                
        #         inGeneKey = "/nfshome/k.bankole/mnt/exasmb.rc.ufl.edu-blue/mcintyre/share/references/dsan_Prin_1.1/dsan11_2_dsan1_gene_key.csv"
        #         inRefAnno = "/nfshome/k.bankole/mnt/exasmb.rc.ufl.edu-blue/mcintyre/share/references/dsan_Prin_1.1/GCF_016746245.2/genomic_gtf_no_gene_features.gtf"
        #         inUJCIndex = "/nfshome/k.bankole/mnt/exasmb.rc.ufl.edu-blue/mcintyre/share/references/dsan_Prin_1.1/dsan11_2_dsan1_ujc_xscript_index.csv"
        #         inUnmapGTF = "/nfshome/k.bankole/mnt/exasmb.rc.ufl.edu-blue/mcintyre/share/references/dsan_Prin_1.1/dsan11_2_dsan1_unmapped.gtf"

        # elif species == "ser":
                
        #         inGeneKey = "/nfshome/k.bankole/mnt/exasmb.rc.ufl.edu-blue/mcintyre/share/references/dser1.1/dser11_2_dser1_gene_key.csv"
        #         inRefAnno = "/nfshome/k.bankole/mnt/exasmb.rc.ufl.edu-blue/mcintyre/share/references/dser1.1/GCF_002093755.2/genomic_gtf_no_gene_features.gtf"
        #         inUJCIndex = "/nfshome/k.bankole/mnt/exasmb.rc.ufl.edu-blue/mcintyre/share/references/dser1.1/dser11_2_dser1_ujc_xscript_index.csv"
        #         inUnmapGTF = "/nfshome/k.bankole/mnt/exasmb.rc.ufl.edu-blue/mcintyre/share/references/dser1.1/dser11_2_dser1_unmapped.gtf"

        # elif species == "yak":
                
        #         inGeneKey = "/nfshome/k.bankole/mnt/exasmb.rc.ufl.edu-blue/mcintyre/share/references/dyak_rogers/dyak105_2_dyak1_gene_key.csv"
        #         inRefAnno = "/nfshome/k.bankole/mnt/exasmb.rc.ufl.edu-blue/mcintyre/share/references/dyak_rogers/Final.DyakMerged_tab_corrected.gtf"
        #         inUJCIndex = "/nfshome/k.bankole/mnt/exasmb.rc.ufl.edu-blue/mcintyre/share/references/dyak_rogers/dyak105_2_dyak1_ujc_xscript_index.csv"
        #         inUnmapGTF = "/nfshome/k.bankole/mnt/exasmb.rc.ufl.edu-blue/mcintyre/share/references/dyak_rogers/dyak105_2_dyak1_unmapped.gtf"
                
        # elif species == "sim":
                
        #         inGeneKey = "/nfshome/k.bankole/mnt/exasmb.rc.ufl.edu-blue/mcintyre/share/references/dsim_WXD1_ASM438218v1/dsimWXD_2_dsimW_gene_key.csv"
        #         inRefAnno = "/nfshome/k.bankole/mnt/exasmb.rc.ufl.edu-blue/mcintyre/share/sex_specific_splicing/dsWXD_union/dsWXD_union_UJC.gtf"
        #         inUJCIndex = "/nfshome/k.bankole/mnt/exasmb.rc.ufl.edu-blue/mcintyre/share/references/dsim_WXD1_ASM438218v1/dsimWXD_2_dsimW_ujc_xscript_index.csv"
        #         inUnmapGTF = "/nfshome/k.bankole/mnt/exasmb.rc.ufl.edu-blue/mcintyre/share/references/dsim_WXD1_ASM438218v1/dsimWXD_2_dsimW_unmapped.gtf"
        
        # Read in the two GTFs
        refGTFDf = trand.io.read_exon_data_from_file(inRefAnno)
        
        try:
                unmapGTFDF = trand.io.read_exon_data_from_file(inUnmapGTF)
        except pd.errors.EmptyDataError:
                print ("Unmapped GTF is empty. Continuing...")
                unmapGTFDF = pd.DataFrame(columns=['transcript_id','gene_id'])
        
        print ("There are " + str(refGTFDf['transcript_id'].nunique()) + " unique transcripts in the reference.")
        print ("There are " + str(unmapGTFDF['transcript_id'].nunique()) + " unmapped transcripts when self-mapping.")
        
        print ("-> " + str(refGTFDf['transcript_id'].nunique() - unmapGTFDF['transcript_id'].nunique()) + " successfully mapped back onto the genome.")
                
        totalMappedTr = (refGTFDf['transcript_id'].nunique() - unmapGTFDF['transcript_id'].nunique())
        
        # Flag Unmapped Transcripts in Reference
        unmapMergeDf = pd.merge(refGTFDf, unmapGTFDF, how='outer', indicator='merge_check')
        
        # Merge Checks that the unmapped is a subset of a reference
        if len(unmapMergeDf) != len(refGTFDf):
                print ("Somehow, there are transcripts in the unmapped GTF that aren't in the reference...")
                # exit()
                
        elif 'right_only' in unmapMergeDf['merge_check']:
                print ("Somehow, there are transcripts in the unmapped GTF that aren't in the reference...")
                # exit()

        # Flag unmapped transcripts
        unmapMergeDf['flagUnmap'] = unmapMergeDf['merge_check'].apply(lambda x: 1 if x == 'both' else 0)
        
        # Convert to simple transcriptID/geneID columns only, with unmap flag
        refDf = unmapMergeDf.drop(columns=['merge_check']).copy(deep=True)        
        refDf = refDf.groupby(['gene_id','transcript_id']).agg("first").reset_index()[['gene_id','transcript_id','flagUnmap']]
        refXscript2Gene = refDf.rename(columns={'gene_id':'ref_geneID','transcript_id':'ref_transcriptID'})

        # Read in gene key from correct_gffcompare output. Verified that unique transcripts in this file matches unique transcripts in corrected_associated...
        keyDf = pd.read_csv(inGeneKey, low_memory=False)
        print ("There are " + str(keyDf['transcript_id'].nunique()) + " unique reference transcripts in the minimap2 gtf/output from correct_gffcompare.")
        
        # Convert to simple transcriptID/geneID columns only
        xscript2AsgnGn = keyDf[['transcript_id','output_gene_id']].copy(deep=True)
        xscript2AsgnGn.columns = ['ref_transcriptID','asgn_geneID']

        refAsgnMergeDf = pd.merge(refXscript2Gene, xscript2AsgnGn, on=['ref_transcriptID'], how='outer',indicator='merge_check')
        
        def checkAsgnGID(mergeCheck, unMap, asgnGID):
                if mergeCheck == 'both':
                        return asgnGID
                elif mergeCheck == 'left_only' and unMap == 1: 
                        return "TR_DID_NOT_MAP"
                elif mergeCheck == 'left_only' and unMap != 1:
                        return "TR_NOT_IN_FASTA"
                elif mergeCheck == 'right_only':
                        return "TR_NOT_IN_REF"
                else:
                        return np.nan
        
        refAsgnMergeDf['asgn_geneID'] =  refAsgnMergeDf.apply(lambda x: checkAsgnGID(x['merge_check'],x['flagUnmap'],x['asgn_geneID']), axis=1)
        
        # Check for transcripts that aren't in the reference's EXON ROWS
        # (if a transcript is only in the corrected GTF, the above situation is the only way i've seen this is possible)
        # When reading the GTF, we only read the rows that have an exon feature. If there is a transcript in the FASTA that
        # does not have any exon rows, it will not be in the list of reference transcripts here in the script.
        
        # so far this only happens in mel650....
        noExonRow = refAsgnMergeDf[(refAsgnMergeDf['asgn_geneID'] == 'TR_NOT_IN_REF')]['ref_transcriptID'].unique()
        if noExonRow.size > 0:
                print("There are {} transcripts do not have exon rows in the reference, but appear in the corrected GTF. Adding to remove list...".format(len(noExonRow)))
                # for tr in noExonRow:
                        # print (tr)
                
                tracker = xscript2AsgnGn[~xscript2AsgnGn['ref_transcriptID'].isin(noExonRow)]
                print ("Value for new number of unique transcripts in the corrected GTF: {}".format(tracker['ref_transcriptID'].nunique()))

                # Remove from output
                refAsgnMergeDf = refAsgnMergeDf[~refAsgnMergeDf['asgn_geneID'].str.contains('TR_NOT_IN_REF')]
                
        # Check for transcripts that aren't in the tr FASTA but are in the reference 
        #(if a transcript is only in the reference but did not map, the above situation is the only way i can think that this is possible)
        notInFASTA = refAsgnMergeDf[refAsgnMergeDf['asgn_geneID'] == 'TR_NOT_IN_FASTA']['ref_transcriptID'].unique() 
        if notInFASTA.size > 0:
                print("There are {} transcripts that not unmapped, but only appear in the reference annotation.".format(len(notInFASTA)))
                
                # for tr in notInFASTA['ref_transcriptID'].unique():
                #         print (tr)
        
                # Remove from output
                refAsgnMergeDf = refAsgnMergeDf[~refAsgnMergeDf['asgn_geneID'].str.contains('TR_NOT_IN_FASTA')]        
                print ("Value for new number of unique reference transcripts: {}".format(refAsgnMergeDf['ref_transcriptID'].nunique()))


        # Verification
        if notInFASTA.size > 0:
                print ("Verification: ")
                print("Number of unique transcripts in reference (excluding those not in the FASTA): {}".format(refGTFDf['transcript_id'].nunique() - len(notInFASTA)))
                print("Number of unique transcripts in output: {}".format(refAsgnMergeDf['ref_transcriptID'].nunique()))
                
                if refGTFDf['transcript_id'].nunique() - len(notInFASTA) == refAsgnMergeDf['ref_transcriptID'].nunique():
                        print("Total number of unique transcripts in output is accurate!")
                else:
                        print("Verification error!")
        else:
                print("Number of unique transcripts in reference (excluding those not in the FASTA): {}".format(refGTFDf['transcript_id'].nunique()))
                print("Number of unique transcripts in output: {}".format(refAsgnMergeDf['ref_transcriptID'].nunique()))
                
                if refGTFDf['transcript_id'].nunique() == refAsgnMergeDf['ref_transcriptID'].nunique():
                        print("Total number of unique transcripts in output is accurate!")
                else:
                        print("Verification error!")
                
        refAsgnMergeDf['flagDiffGeneID'] = refAsgnMergeDf.apply(lambda x: 1 if x['ref_geneID'] != x['asgn_geneID'] else 0,axis=1)
        
        
        #Output.
        outDf = refAsgnMergeDf.drop(columns='merge_check')
        outDf = outDf[['ref_transcriptID','ref_geneID','asgn_geneID','flagDiffGeneID','flagUnmap']]
        
        outFile = "{}/{}_gffcompare_correction_summary.csv".format(outDir, filePrefix)
        outDf.to_csv(outFile,index=False)
        
        
        # Output a list of transcripts to remove from the corrected GTF!
        removeFile = '{}/{}_remove_these_transcripts_from_corrected_gtf.txt'.format(outDir, filePrefix)
        with open(removeFile, 'w+') as f:
                removeList = noExonRow  
                f.write("\n".join(removeList))     

if __name__ == '__main__':
        # # Parse command line arguments
        global args
        args = getOptions()
        main()
