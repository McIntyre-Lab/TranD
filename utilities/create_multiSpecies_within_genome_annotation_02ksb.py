#!/usr/bin/env python

import argparse
import glob
import os
import pandas as pd
import trand.io


def getOptions():
    # Parse command line arguments
    parser = argparse.ArgumentParser(description=
                                     "Input genomeName as seen in mapping/ujc output, "
                                     "path to UJC GTF files all mapped to said genome in one argument, and"
                                     "outfile for the GTF (path and name). Outputs a file only contains unique jxnHashes "
                                     "from all species within one genome.")

    # Input data
    parser.add_argument("-gn",
                        "--genome-name",
                        dest="genomeName",
                        required=True,
                        help="Name of the coordinates that all the UJC indexes output are mapped to."
                        "Must match the file name (ex: example_2_{genomeName}_ujc.gtf")

    parser.add_argument("-ig",
                        "--input-gtfs",
                        dest="inCat",
                        nargs="+",
                        required=True,
                        help="All paths to location of UJC GTFs mapped to same genome.")

    # Output data
    parser.add_argument("-o",
                        "--output-file",
                        dest="outGTF",
                        required=True,
                        help="Path and name of GTF to output to (ex: "
                        "/path/to/output/fiveSpecies_2_dmel6_ujc.gtf")

    args = parser.parse_args()
    return args


def main():
    # Parse command line arguments

    genomeName = "dmel6"
    inCatLst = ["/nfshome/k.bankole/mnt/exasmb.rc.ufl.edu-blue/mcintyre/share/sex_specific_splicing/annotation_id_ujc_output/dmel650_2_dmel6_ujc.gtf",
             "/nfshome/k.bankole/mnt/exasmb.rc.ufl.edu-blue/mcintyre/share/sex_specific_splicing/mapped_ujc_id_ujc_output/dsim202_2_dsim2_ujc_2_dmel6_noGeneID_ujc.gtf",
             "/nfshome/k.bankole/mnt/exasmb.rc.ufl.edu-blue/mcintyre/share/sex_specific_splicing/mapped_ujc_id_ujc_output/dsan11_2_dsan1_ujc_2_dmel6_noGeneID_ujc.gtf",
             "/nfshome/k.bankole/mnt/exasmb.rc.ufl.edu-blue/mcintyre/share/sex_specific_splicing/mapped_ujc_id_ujc_output/dsimWXD_2_dsim2_ujc_2_dmel6_noGeneID_ujc.gtf",
             "/nfshome/k.bankole/mnt/exasmb.rc.ufl.edu-blue/mcintyre/share/sex_specific_splicing/mapped_ujc_id_ujc_output/dyak21_2_dyak2_ujc_2_dmel6_noGeneID_ujc.gtf",
             "/nfshome/k.bankole/mnt/exasmb.rc.ufl.edu-blue/mcintyre/share/sex_specific_splicing/mapped_ujc_id_ujc_output/dser11_2_dser1_ujc_2_dmel6_noGeneID_ujc.gtf"]

    outGTF = "/nfshome/k.bankole/mnt/exasmb.rc.ufl.edu-blue/mcintyre/share/sex_specific_splicing/fiveSpecies_annotations/fiveSpecies_2_dmel6_ujc_roz.gtf"

    # genomeName = args.genomeName
    # inCatLst = args.inCat
    # outGTF = args.outGTF


    # TODO: keep track of number of exons for each transcript to check    
        
    gtfDct = dict()
    for file in inCatLst:
        #do things
        fileName = os.path.basename(file)
        print("Reading file: {}".format(fileName))
    
        annoName = fileName.split('_2')[0].rstrip('_')
        df = trand.io.read_exon_data_from_file(file)
        gtfDct[annoName] = df
    
    
    allHashSet = set()
    masterDf = pd.DataFrame()
    masterDf = list(gtfDct.values())[0]
    gtf = list(gtfDct.values())[1]
    
    testDf = masterDf[masterDf['transcript_id'] == "000ab5cb90829f133fd4e78f093a2d795bc47c83a7aaef93743560cf7c37f949"]
    test2Df = gtf[gtf['transcript_id'] == "000ab5cb90829f133fd4e78f093a2d795bc47c83a7aaef93743560cf7c37f949"]
    
    # group by transcript and merge............
    merge = pd.merge(testDf, test2Df, on=['seqname','start','end','strand','gene_id','transcript_id'], how='outer', indicator='merge_check').drop_duplicates()
    
    # TODO: check that theres only two in left and right only...
    allGood = merge[merge['merge_check'] == "both"].drop('merge_check', axis=1)
    
    leftOnly = merge[merge['merge_check'] == "left_only"].drop('merge_check', axis=1)
    rightOnly = merge[merge['merge_check'] == "left_only"].drop('merge_check', axis=1)

    leftOnly[]

    leftOnly = leftOnly.agg({
        'start':max,
        'end':max})


    masterDf.groupby('transcript_id')
    for row in masterDf.to_dict('records'):
        jHash = row['transcript_id']
        
        if jHash in ['gtf']
        
        
        

    
    testDf = masterDf[masterDf['transcript_id'] == "000ab5cb90829f133fd4e78f093a2d795bc47c83a7aaef93743560cf7c37f949"]
    test2Df = gtf[gtf['transcript_id'] == "000ab5cb90829f133fd4e78f093a2d795bc47c83a7aaef93743560cf7c37f949"]
    


    
    for gtf in gtfDct.values():
        if masterDf.empty:
            masterDf = gtf
        else:
            
            # TODO: verify num in both, num left only, num right only is accurate
            
            # TODO: verify that this works
            merge = pd.merge(masterDf, gtf, on=['seqname','strand','start','end','gene_id','transcript_id'], how='outer', indicator='merge_check').drop_duplicates()
            
            merge['start'] = merge[['start_x','start_y']].max(axis=1)
            merge.drop(['start_x','start_y'], inplace=True, axis=1)

            merge['end'] = merge[['end_x','end_y']].max(axis=1)
            merge.drop(['end_x','end_y'], inplace=True, axis=1)
            
            masterDf = merge.drop_duplicates().drop('merge_check', axis = 1)    

    masterDf['transcript_id'].nunique()
    
    masterDf.to_csv(index=False)
    
    print("Writing GTF...")
    if os.path.isfile(outGTF):
        os.remove(outGTF)

    trand.io.write_gtf(data=masterDf, out_fhs={"gtf": "/nfshome/k.bankole/mclab/SHARE/McIntyre_Lab/sex_specific_splicing/one_gene_her_trand_test/test.gtf"}, fh_name="gtf")
    
    # Create dct {ujc gtf name:GTF DF}
    # indexDfDct = dict()
    # allHashSet = set()

    # for file in glob.glob(inFolder + '/*_2_{}_noGeneID_ujc.gtf'.format(genomeName)) + glob.glob(inFolder + '/*_2_{}_ujc_roz.gtf'.format(genomeName)):

    #     fileName = os.path.basename(file)
    #     print("Reading file: {}".format(fileName))

    #     annoName = fileName.split('2')[0].rstrip('_')
    #     inDf = trand.io.read_exon_data_from_file(file)
    #     # print(len(allHashSet))

    #     indexDfDct[annoName] = inDf
    #     allHashSet.update(set(inDf['transcript_id'].unique().tolist()))
    #     # print(len(allHashSet))

    # print("Total Number of Unique UJCs across all annotations: {}".format(
    #     len(allHashSet)))

    # # test = set()
    # # for annoName, df in list(indexDfDct.items()):
    # #     # print(annoName)
    # #     test.update(set(df['transcript_id'].unique().tolist()))

    # # len(set(test))

    # annoName = list(indexDfDct.keys())[0]
    # concatDf = indexDfDct[annoName].copy(deep=True)

    # for annoName, df in list(indexDfDct.items())[1:]:
    #     concatDf = pd.concat(
    #         [concatDf, df]).drop_duplicates().reset_index(drop=True)

    # print("Total Number of Unique UJCs in output: {}".format(
    #     concatDf['transcript_id'].nunique()))

    # # Could use some extra verification but... seems like this is good?
    # # You could not have been more wrong^^^^^^^^

    # print("Writing GTF...")
    # if os.path.isfile(outGTF):
    #     os.remove(outGTF)

    # trand.io.write_gtf(data=concatDf, out_fhs={"gtf": outGTF}, fh_name="gtf")

    # print("Complete!")


if __name__ == '__main__':
    global args
    args = getOptions()
    main()
