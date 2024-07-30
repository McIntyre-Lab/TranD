# #!/usr/bin/env python

# import argparse
# import trand.io
# import time
# import pandas as pd


# def getOptions():
#     # Parse command line arguments
#     parser = argparse.ArgumentParser(description="")

#     # Input data
#     parser.add_argument("-", "--", dest="", required=True, help="")

#     # Output data
#     parser.add_argument("-", "--", dest="", required=True, help="")

#     args = parser.parse_args()
#     return args


# def read_exon_data_from_file(infile):
#     print("Reading GTF...")

#     omegatic = time.perf_counter()

#     all_gtf_columns = ['seqname', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame',
#                        'attributes', 'comments']

#     drop_columns = ['source', 'feature', 'score', 'frame', 'comments']
#     drop_columns = ['score', 'frame', 'comments']

#     data = pd.read_csv(infile, sep='\t', comment='#',
#                        header=None, low_memory=False)
#     file_cols = data.columns

#     if len(file_cols) < len(all_gtf_columns):
#         gtf_cols = all_gtf_columns[:len(file_cols)]
#     data.columns = gtf_cols
#     drop_cols = [x for x in drop_columns if x in gtf_cols]

#     data = data[data['feature'].isin(['gene', 'exon'])]
#     data = data.drop(labels=drop_cols, axis=1)

#     data['seqname'] = data['seqname'].astype(str)
#     data['start'] = data['start'].astype(int)
#     data['end'] = data['end'].astype(int)

#     data.reset_index(drop=True, inplace=True)

#     seqnameLst = []
#     startLst = []
#     endLst = []
#     strandLst = []
#     geneIDLst = []
#     xscriptIDLst = []
#     featureLst = []

#     for row in data.to_dict('records'):
#         feature = row['feature']
#         rawAttr = row['attributes']
#         attrLst = [x.strip() for x in rawAttr.strip().split(';')]
#         gnTrAttr = [x for x in attrLst if (
#             'transcript_id' in x and 'orig_transcript_id' not in x) or 'gene_id' in x]
#         gene_id, transcript_id = None, None

#         for item in gnTrAttr:
#             if 'gene_id' in item:
#                 gene_id = item.split('gene_id')[1].strip().strip('\"')
#             elif 'transcript_id' in item and 'orig_transcript_id' not in item:
#                 transcript_id = item.split(
#                     'transcript_id')[-1].strip().strip('\"')

#         if not gene_id and row['feature'] == 'gene':
#             print("gene_id not found in '{}'", row)
#             gene_id = None

#         if not transcript_id and row['feature'] != 'gene':
#             print("transcript_id not found in '{}'", row)
#             transcript_id = None

#         seqnameLst.append(row['seqname'])
#         startLst.append(row['start'])
#         endLst.append(row['end'])
#         strandLst.append(row['strand'])

#         geneIDLst.append(gene_id)
#         xscriptIDLst.append(transcript_id)
#         featureLst.append(feature)

#     newData = pd.DataFrame(
#         {
#             'feature': featureLst,
#             'seqname': seqnameLst,
#             'start': startLst,
#             'end': endLst,
#             'strand': strandLst,
#             'gene_id': geneIDLst,
#             'transcript_id': xscriptIDLst
#         })

#     missing_value_num = newData.isnull().sum().sum()
#     if missing_value_num > 0:
#         raise Exception("whoopsie :3")
#     else:
#         print("No missing values in data")

#     gene_id_missing_value_num = newData['gene_id'].isnull().sum()

#     transcript_id_missing_value_num = newData['transcript_id'].isnull().sum()

#     if gene_id_missing_value_num > 0:
#         raise Exception("whoopsie :3")
#     if transcript_id_missing_value_num > 0:
#         print("Missing transcript_id value number: {}",
#               transcript_id_missing_value_num)

#     newData['start'] = pd.to_numeric(newData['start'], downcast="unsigned")
#     newData['end'] = pd.to_numeric(newData['end'], downcast="unsigned")

#     toc = time.perf_counter()

#     print(f"GTF Read complete,  took {toc-omegatic:0.4f} seconds.")
#     return newData


# refGTF = "/nfshome/k.bankole/mnt/exasmb.rc.ufl.edu-blue/mcintyre/share/references/dsan_Prin_1.1/dsan11_protein_coding_ref.gtf"

# dataGTF = "/nfshome/k.bankole/Desktop/gffcompare_replacement.gtf"


# refDf = read_exon_data_from_file(refGTF)
# refDf = refDf[refDf['feature'] == 'gene']

# dataDf = read_exon_data_from_file(dataGTF)

# dataDf['transcript_id'] = dataDf['transcript_id'].astype(str)
# dataDf['gene_id'] = dataDf['gene_id'].astype(str)

# xscriptDf = dataDf.groupby('transcript_id').agg({
#     'seqname': set,
#     'start': min,
#     'end': max,
#     'strand': set
# }).reset_index()

# singleChrXscript = xscriptDf['seqname'].apply(lambda x: len(x) == 1)
# if not singleChrXscript.all():
#     raise Exception(
#         "There are ERPs belonging to more than one chromosome. Quitting.")
# else:
#     xscriptDf['seqname'] = xscriptDf['seqname'].apply(
#         lambda x: list(x)[0])

# singleStrandXscript = xscriptDf['strand'].apply(lambda x: len(x) == 1)
# if not singleStrandXscript.all():
#     raise Exception(
#         "There are ERPs belonging to more than one strand. Quitting.")
# else:
#     xscriptDf['strand'] = xscriptDf['strand'].apply(
#         lambda x: list(x)[0])


# refDf = refDf.sort_values(['seqname', 'start'])
# refDf = refDf.drop('transcript_id', axis=1)
# xscriptDf = xscriptDf.sort_values(['seqname', 'start'])


# mergedDf = pd.merge(xscriptDf, refDf, on=[
#                     'seqname', 'strand'], suffixes=('_xscript', '_ref'))


# def check_overlap(row):
#     return max(row['start_xscript'], row['start_ref']) < min(row['end_xscript'], row['end_ref'])


# # Apply the overlap check
# mergedDf['overlap'] = mergedDf.apply(check_overlap, axis=1)

# matched = mergedDf[mergedDf['overlap']]
# unmatched = mergedDf[~mergedDf['overlap']]

# omegatic = time.perf_counter()

# geneDct = dict()
# unmatchLst = []
# for row in xscriptDf.to_dict('records'):

#     chromosome = row['seqname']
#     strand = row['strand']
#     xscript = row['transcript_id']
#     xscriptStart = row['start']
#     xscriptEnd = row['end']

#     for gnRow in refDf.to_dict('records'):
#         if chromosome != gnRow['seqname'] or strand != gnRow['strand']:
#             continue
#         else:
#             if max(xscriptStart, gnRow['start']) < min(xscriptEnd, gnRow['end']):
#                 # print(row)
#                 # print(erID)
#                 # print(erInfo)

#                 if gnRow['gene_id'] in geneDct.keys():
#                     geneDct[gnRow['gene_id']].append(xscript)
#                 else:
#                     geneDct[gnRow['gene_id']] = [xscript]

#                 break
#     else:
#         unmatchLst.append(xscript)


# toc = time.perf_counter()

# print(f"Gene overlappign complete, took {toc-omegatic:0.4f} seconds.")


# # if __name__ == '__main__':
# #     # Parse command line arguments
# #     global args
# #     args = getOptions()
# #     main()
