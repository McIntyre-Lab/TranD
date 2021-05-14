#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
# Copyright © 2021 Oleksandr Moskalenko <om@rc.ufl.edu>
# Distributed under terms of the MIT license.
"""
Event Analysis 2.0 is an evolution of the concepts in the original EA from
https://github.com/McIntyre-Lab/events, but simplified for the task of comparing data produced by
sequence alignment and annotation tools.

This program takes a single GTF file for all pairwise transcript comparisons within each gene or two
GTF files for transcript comparisions between genes in both datasets, reads in exonic fragment data,
and performs event analysis to determine transcript distance measures for the genes and transcripts
in the data.

Raw gtf data:
seqname source  feature  start end  score strand frame attributes
2L      FlyBase 5UTR     7529  7679 .     +      .     gene_symbol "CG11023"; transcript_id "FBtr...
2L      FlyBase exon     7529  8116 .     +      .     gene_symbol "CG11023"; transcript_id "FBtr...

Exon Fragment Data for Event Analysis:
seqname start end   strand gene_id      transcript_id
2L      7529  8116  +      FBgn0031208  FBtr0300689
2L      8193  9484  +      FBgn0031208  FBtr0300689
2L      9839  11344 -      FBgn0002121  FBtr0078169

"""

import argparse
import itertools
import logging
import sys
import pandas as pd
from collections import namedtuple
from loguru import logger
from numpy import nan
from pathlib import Path
from pybedtools import BedTool

# Import transcript distance functions and variables
import transcript_distance as TD

# DATA OUTPUT CONFIGURATION
common_outfiles = {'ea_fh': 'event_analysis.csv', 'jc_fh': 'junction_catalog.csv'}

ea_df_cols = ['gene_id', 'transcript_1', 'transcript_2', 'transcript_id', 'ef_id', 'ef_chr',
              'ef_start', 'ef_end', 'ef_strand', 'ef_ir_flag', 'er_id', 'er_chr', 'er_start',
              'er_end', 'er_strand']
jct_df_cols = ['gene_id', 'transcript_id', 'coords']

ea_df_gene_cols = ['gene_id', 'transcript_1', 'transcript_2', 'transcript_id', 'ef_id', 'ef_chr',
              'ef_start', 'ef_end', 'ef_strand', 'ef_ir_flag', 'er_id', 'er_chr', 'er_start',
              'er_end', 'er_strand']

pairwise_outfiles = {'td_fh':'pairwise_transcript_distance.csv'}
# Later when TD has been added
# common_outfiles = {'ea_er_fh': 'gene_ea_exonic_regions.csv', 'ea_ef_fh':
#                    'gene_ea_exonic_fragments.csv', 'ea_fh': 'pairwise_ea.csv', 'td_fh':
#                    'transcript_distances.csv'}

two_gtfs_outfiles = {'gtf1_fh': 'gtf1_only.gtf', 'gtf2_fh': 'gtf2_only.gtf'}
# , 'gtf_names_fh': 'gtf_metadata.csv'}


def parse_args(print_help=False):
    """Parse command-line arguments"""
    class MyParser(argparse.ArgumentParser):
        """Subclass ArgumentParser for better help printing"""
        def error(self, message):
            sys.stderr.write("error: %s\n" % message)
            self.print_help()
            sys.exit(2)
    parser = MyParser(
        description="Perform Event Analysis on annotated NGS read data."
    )
    parser.add_argument(
        dest="infiles",
        metavar="input_file",
        type=str,
        nargs='+',
        help="One or two input GTF (or GFFv2) file(s).",
    )
    parser.add_argument(
        "-o",
        "--outdir",
        action="store",
        type=str,
        required=False,
        help="Output directory. Default:current directory.",
    )
    parser.add_argument(
        "-l",
        "--logfile",
        dest="log_file",
        action="store",
        type=str,
        default=None,
        required=False,
        help="Log file name for logging processing events to file.",
    )
    parser.add_argument(
        "-e", "--ea",
        dest='ea_mode',
        type=str,
        choices=['pairwise', 'gene'],
        default='pairwise',
        help="Event analysis based on a pair of transcripts for TD or all gene isoforms without TD",
    )
    parser.add_argument(
        "-f", "--force",
        action="store_true",
        help="Force overwrite existing output directory and files within.",
    )
    parser.add_argument(
        "-v", "--verbose",
        action="store_true",
        help="Verbose output",
    )
    parser.add_argument(
        "-d", "--debug",
        action="store_true",
        help=argparse.SUPPRESS
    )
    if print_help:
        parser.print_help()
        sys.exit(0)
    args = parser.parse_args()
    if len(args.infiles) > 2:
        print("\nToo many input files - pass one or two GTF/GFF files as input.\n")
        parser.print_help()
        sys.exit(2)
    if args.ea_mode == 'gene':
        if len(args.infiles) > 1:
            logger.warning("EA 'gene' mode is ignored for two GTF files - only pairwise is done.")
    return args


def handle_outdir(args):
    if not args.outdir:
        outdir = Path.cwd()
    else:
        outdir = Path(args.outdir)
    logger.debug("Output directory: {}", str(outdir))
    if outdir.exists():
        if not args.force:
            exit("Not overwriting existing output directory without -f|--force. Exiting.")
    outdir.mkdir(parents=True, exist_ok=True)
    return(str(outdir))


def setup_logging(debug, verbose, logfile):
    """Set the correct logging level and sinks."""
    logger.remove()
    if verbose:
        level = logging.INFO
    else:
        level = logging.WARN
    if debug:
        level = logging.DEBUG
        logger.add(sys.stderr, level=level)
        logger.debug("Debugging output enabled")
        if logfile:
            logger.debug("Logging to {}", logfile)
    else:
        logger.add(sys.stderr, level=level)
        if logfile:
            logger.info("Logging to {}", logfile)
    if logfile:
        logger.add(logfile, level=level)
    logger.debug("Logging level set to : {}", level)


def open_output_files(outdir, outfiles):
    """Open all necessary output files and return filehandles"""
    out_fhs = {}
    try:
        for outfile in outfiles:
            fh = Path(outdir) / outfiles[outfile]
            open(fh, 'w')
            out_fhs[outfile] = fh
    except SystemError as e:
        logger.error("Cannot open one or more output file for writing: {}", e)
    return out_fhs


def read_exon_data_from_file(infile):
    """Create a pandas dataframe with exon records from a gtf file"""
    all_gtf_columns = ['seqname', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame',
                       'attributes', 'comments']
    drop_columns = ['source', 'feature', 'score', 'frame', 'comments']
    data = pd.read_csv(infile, sep='\t', comment='#', header=None, low_memory=False)
    file_cols = data.columns
    if len(file_cols) < len(all_gtf_columns):
        gtf_cols = all_gtf_columns[:len(file_cols)]
    data.columns = gtf_cols
    drop_cols = [x for x in drop_columns if x in gtf_cols]
    logger.debug("Raw data rows: {}", data.shape[0])
    data = data[data['feature'] == 'exon']
    data = data.drop(labels=drop_cols, axis=1)
    data.reset_index(drop=True, inplace=True)
    for i in data.index:
        raw_attrs = data.at[i, "attributes"]
        attr_list = [x.strip() for x in raw_attrs.strip().split(';')]
        g_t_attrs = [x for x in attr_list if 'transcript_id' in x or 'gene_id' in x]
        gene_id, transcript_id = nan, nan
        for item in g_t_attrs:
            if 'gene_id' in item:
                gene_id = item.split('gene_id')[1].strip().strip('\"')
            elif 'transcript_id' in item:
                transcript_id = item.split('transcript_id')[-1].strip().strip('\"')
        if not gene_id:
            logger.error("gene_id not found in '{}'", data[i])
        if not transcript_id:
            logger.error("transcript_id not found in '{}'", data[i])
        # logger.debug("Gene: {}, Transcript: {}", gene_id, transcript_id)
        data.at[i, "gene_id"] = gene_id
        data.at[i, "transcript_id"] = transcript_id
    data = data.drop(labels='attributes', axis=1)
    logger.debug("Exon data rows: {}", data.shape[0])
    missing_value_num = data.isnull().sum().sum()
    if missing_value_num > 0:
        logger.warning("Total number of missing values: {}", missing_value_num)
    else:
        logger.info("No missing values in data")
    gene_id_missing_value_num = data['gene_id'].isnull().sum()
    transcript_id_missing_value_num = data['transcript_id'].isnull().sum()
    if gene_id_missing_value_num > 0:
        logger.warning("Missing gene_id value number: {}", missing_value_num)
    if transcript_id_missing_value_num > 0:
        logger.warning("Missing transcript_id value number: {}", missing_value_num)
    data['start'] = pd.to_numeric(data['start'], downcast="unsigned")
    data['end'] = pd.to_numeric(data['end'], downcast="unsigned")
    # str_columns = ['seqname', 'strand', 'gene_id', 'transcript_id']
    # for col in str_columns:
    #     data[col] = data[col].astype('|S')
    # logger.debug(data.memory_usage(deep=True))
    # logger.debug("Data types: {}", data.dtypes)
    return data


def verify_same_reference(tx_names, data):
    """Verify that exons are on the same contig/chrom"""
    locs = data['seqname'].unique()
    if len(locs) > 1:
        raise ValueError(f"Multiple contig/chromn locations found for {tx_names}, skipping.")


def verify_same_strand(tx_names, data):
    """Verify exons are on the same strand"""
    strands = data['strand'].unique()
    if len(strands) > 1:
        raise ValueError(f"Multiple strands found in {tx_names}, skipping.")


def write_td_data(data, out_fhs):
    """Write results of a transcript distance comparison to output files."""
    logger.debug("Writing TD data")


def do_td(data):
    """Perform Transcript Distance analysis on a pair of transcripts"""
    raise NotImplementedError
    # tx_names = data['transcript_id'].unique()
    # logger.debug("TD on {}", tx_names)
    td_data = ""
    return td_data


def write_output(data, out_fhs, fh_name):
    """Write results of event analysis to output files."""
    data.to_csv(out_fhs[fh_name], mode='a', header=False, index=False)


def prep_bed_for_ea(data):
    """Event analysis"""
    exons = {}
    gene_id = data['gene_id'].unique()[0]
    tx_names = list(data['transcript_id'].unique())
    try:
        verify_same_strand(tx_names, data)
        verify_same_reference(tx_names, data)
    except ValueError:
        raise
    exons = {'gene_id': gene_id}
    exons['transcript_list'] = tx_names
    for tx in tx_names:
        exon_id = 1
        tx_data = data[data['transcript_id'] == tx]
        tx_bed = []
        for index, row in tx_data.iterrows():
            # Score is optional, use zero to keep the strand data in the next column
            tx_bed.append((row.seqname, str(int(row.start)-1), str(int(row.end)),
                           f"{row.transcript_id}_exon_{str(exon_id)}", '0', row.strand))
            exon_id += 1
        exons[tx] = tx_bed
    return exons


def format_ea_identical(gene, tx_names, data):
    """Format exons as EFs for identical transcripts"""
    # Start ef_id numbering from 1
    # ef_id = 1
    for feature in data:
        pass
    output = data
    return(output)


def create_junction_catalog(gene, tx, tx_data):
    """Create a junction catalog for a transcript"""
    junctions = []
    id = 0
    for e in tx_data:
        id += 1
        if id == 1:
            left_end = int(e.end) - 10
            continue
        right_start = int(e.start + 10)
        jct_coords = f"{e.chrom}:{left_end}:{right_start}:{e.strand}"
        junctions.append([gene, tx, jct_coords])
        left_end = int(e.end) - 10
    return junctions

def create_junction_catalog_str(gene, tx, tx_data):
    """Create a junction catalog for a transcript"""
    junctions = []
    id = 0
    tx_data_df = pd.DataFrame(tx_data,columns=['chrom','start','end','name','score','strand'])
    tx_data_df['start'] = tx_data_df['start'].astype(int)
    tx_data_df['end'] = tx_data_df['end'].astype(int)
    tx_data_df = tx_data_df.sort_values('start')
    for index,e in tx_data_df.iterrows():
        id += 1
        if id == 1:
            left_end = int(e.end) - 10
            continue
        right_start = int(e.start + 10)
        jct_coords = f"{e.chrom}:{left_end}:{right_start}:{e.strand}"
        junctions.append([gene, tx, jct_coords])
        left_end = int(e.end) - 10
    return junctions

def get_intron_retention_efs(ers_bed, efs_bed, common_efs):
    """
    Produce a list of exon fragments that participate in intron retention events.
    In essence, retained introns are ER-internal transcript-specific exonic fragments.
    So, they cannot be on the ER borders and cannot be shared between two transcripts.
    For a tiny speedup check that we have more than two EFs in an ER.
    """
    ir_efs = []
    for er in ers_bed:
        er_bed_str = f"{er.chrom}\t{er.start}\t{er.end}\t{er.name}\t{er.score}\t{er.name}"
        er_bed = BedTool(er_bed_str, from_string=True)
        er_efs = efs_bed.intersect(er_bed, wb=True)
        if len(er_efs) >= 3:
            for ef in er_efs:
                # Intron cannot be on the ER border by definition
                if ef.start != er.start and ef.end != er.end:
                    if ef.name not in common_efs:
                        ir_efs.append(ef.name)
    return ir_efs


def do_ea_gene(data):
    """
    Event Analysis on a pair of transcripts
    """
    gene_id = data['gene_id']
    logger.debug("Gene EA on {}", gene_id)
    logger.debug(data)
    tx_names = data['transcript_list']
    tx_data = {}
    for tx_name in tx_names:
        tx_bed_str = data[tx_name]
        logger.debug("Gene EA Raw Tx: {}\n{}", tx_name, tx_bed_str)
        tx_data[tx_name] = BedTool(tx_bed_str).saveas()
    tx0 = tx_names[0]
    junction_data = create_junction_catalog(gene_id, tx0, tx_data[tx0])
    for tx_name in tx_names[1:]:
        junction_data.extend(create_junction_catalog(gene_id, tx_name, tx_data[tx_name]))
    junction_df = pd.DataFrame(junction_data, columns=jct_df_cols)
    logger.debug(tx_data)
    ea_df = pd.DataFrame(columns=ea_df_gene_cols)
    return ea_df, junction_df


def do_ea_pair(data):
    """
    Event Analysis on a pair of transcripts
    """
    gene_id = data['gene_id']
    tx_names = data['transcript_list']
    tx1_name, tx2_name = tx_names[0], tx_names[1]
    tx1_bed_str = data[tx1_name]
    tx2_bed_str = data[tx2_name]
    logger.debug("Comparing {} vs {}", tx1_name, tx2_name)
    junction_data1 = create_junction_catalog_str(gene_id, tx1_name, tx1_bed_str)
    junction_data2 = create_junction_catalog_str(gene_id, tx2_name, tx2_bed_str)
    junction_str1 = "|".join(pd.DataFrame(junction_data1,columns=['gene_id','transcript_id','junction_id'])['junction_id'])
    junction_str2 = "|".join(pd.DataFrame(junction_data2,columns=['gene_id','transcript_id','junction_id'])['junction_id'])
    junction_data = junction_data1 + junction_data2
    # Check if transcripts are full-splice match (all junctions are the same)
    same_start = False
    same_end = False
    is_fsm = False
    if junction_str1 == junction_str2:
        is_fsm = True
        # Check if start coordinates of the transcripts the same
        if pd.DataFrame(tx1_bed_str)[1].min() == pd.DataFrame(tx2_bed_str)[1].min():
            same_start = True
        # Check if end coordinates of transcripts are also the same
        if pd.DataFrame(tx1_bed_str)[2].max() == pd.DataFrame(tx2_bed_str)[2].max():
            same_end = True
        if same_start and same_end:
            # Transcripts are identical
            logger.info("Exons are identical for {}.", " / ".join(tx_names))
            ea_data = format_fsm_pair_ea(tx1_bed_str, tx2_bed_str, tx1_name, tx2_name, gene_id, side_diff="none")
        elif same_start and not same_end:
            # Start coordinates and junctions are the same
            logger.info("Junctions and start coordinates are identical for {}.", " / ".join(tx_names))
            ea_data = format_fsm_pair_ea(tx1_bed_str, tx2_bed_str, tx1_name, tx2_name, gene_id, side_diff="end")
            # !!! Add calculation of distances (will only need to calculate the difference of 3' ends
        elif same_end and not same_start:
            # End coordinates and junctions are the same
            logger.info("Junctions and end coordinates are identical for {}.", " / ".join(tx_names))
            ea_data = format_fsm_pair_ea(tx1_bed_str, tx2_bed_str, tx1_name, tx2_name, gene_id, side_diff="start")
            # !!! Add calculation of distances (will only need to calculate the difference of 5' ends
        else:
            # Junctions are the same and ends are both different
            logger.info("Junctions are identical for {}.", " / ".join(tx_names))
            ea_data = format_fsm_pair_ea(tx1_bed_str, tx2_bed_str, tx1_name, tx2_name, gene_id, side_diff="both")
            # !!! Add calculation of distances (will only need to calculate the difference of 5' and 3' ends)
    else:
        # Junctions are not identical (there is some alternate donor/acceptor/exon)
        tx1_bed = BedTool(tx1_bed_str).saveas()
        tx2_bed = BedTool(tx2_bed_str).saveas()
        ea_data = er_ea_analysis(tx1_bed, tx2_bed, tx1_name, tx2_name, gene_id)
    out_df = pd.DataFrame(ea_data, columns=ea_df_cols)
    junction_df = pd.DataFrame(junction_data, columns=jct_df_cols)
    td_df = pd.DataFrame(TD.calculate_distance(out_df, junction_df, gene_id, tx1_name, tx2_name, fsm=is_fsm)).T
    return out_df, junction_df, td_df


def er_ea_analysis(tx1_bed, tx2_bed, tx1_name, tx2_name, gene_id):
    """Convert ERs (Exonic Regions) into EFs (Exonic Fragments) and analyze all er events"""
    ea_data = []
    strand = list(set([x.strand for x in tx1_bed]))[0]
    # logger.debug("TX1: {}: \n{}", tx1_name, tx1_bed)
    # logger.debug("TX2: {}: \n{}", tx2_name, tx2_bed)
    raw_ers_bed = tx1_bed.cat(tx2_bed, postmerge=True)
    raw_ers_list = [str(x).split() for x in raw_ers_bed]
    er_id = 1
    ers_list = []
    for i in raw_ers_list:
        i.extend([f"{gene_id}:ER{er_id}", '0', strand])
        ers_list.append("\t".join(i))
        er_id += 1
    ers_str = "\n".join(ers_list)
    # logger.debug("ERs for {} and {}: \n{}", tx1_name, tx2_name, ers_str)
    # logger.debug("ERs str:\n{}", ers_str)
    ers_bed = BedTool(ers_str, from_string=True)
    er_data = {}
    er_datum = namedtuple('ER', ['chrom', 'start', 'end', 'strand'])
    for er in ers_bed:
        er_data[er.name] = er_datum(er.chrom, er.start, er.end, er.strand)
    ef_tx1 = tx1_bed.subtract(tx2_bed)
    ef_tx2 = tx2_bed.subtract(tx1_bed)
    ef_both = tx1_bed.intersect(tx2_bed)
    efs_raw = ef_tx1.cat(ef_tx2, postmerge=False).cat(ef_both, postmerge=False).sort()
    # logger.debug("Raw EFs: \n{}", efs_raw)
    efs_er = ers_bed.intersect(efs_raw).sort()
    # logger.debug("EFs in ERs: \n{}", efs_er)
    efs_list = []
    ef_id = 1
    er_name = ''
    for i in efs_er:
        if i.name != er_name:
            ef_id = 1
        er_name = i.name
        ef_name = f"{er_name}:EF{ef_id}"
        ef_str = f"{i.chrom}\t{i.start}\t{i.end}\t{ef_name}\t {i.score}\t{i.strand}"
        efs_list.append(ef_str)
        ef_id += 1
    efs_str = "\n".join(efs_list)
    # logger.debug("Final EFs: \n{}", efs_str)
    efs_bed = BedTool(efs_str, from_string=True)
    common_efs_set = set([x.name for x in efs_bed.intersect(ef_both)])
    common_efs = list(common_efs_set)
    tx1_efs = list(set([x.name for x in efs_bed.intersect(tx1_bed)]).difference(common_efs_set))
    tx2_efs = list(set([x.name for x in efs_bed.intersect(tx2_bed)]).difference(common_efs_set))
    ir_efs = get_intron_retention_efs(ers_bed, efs_bed, common_efs)

    for ef in efs_bed:
        ir_flag = '0'
        ef_name = ef.name
        er_name = ":".join(ef_name.split(':')[:-1])
        if ef_name in ir_efs:
            ir_flag = '1'
        if ef_name in common_efs:
            tx_list = f"{tx1_name}|{tx2_name}"
        elif ef_name in tx1_efs:
            tx_list = f"{tx1_name}"
        elif ef_name in tx2_efs:
            tx_list = f"{tx2_name}"
        else:
            logger.error(f"Cannot find the fragment {ef_name} in any fragment lists, skipping")
        ea_datum = [gene_id, tx1_name, tx2_name, tx_list, ef_name, ef.chrom, ef.start, ef.end,
                    ef.strand, ir_flag, er_name, er_data[er_name].chrom, er_data[er_name].start,
                    er_data[er_name].end, er_data[er_name].strand]
        ea_data.append(ea_datum)
    logger.debug("Final EA data: \n{}", ea_data)
    return ea_data


def format_identical_pair_ea(tx1_bed, tx2_bed, tx1_name, tx2_name, gene_id):
    ea_data = []
    er_id, ef_id = 1, 1
    for i in tx1_bed:
        er_name = f"{gene_id}:ER{er_id}"
        ef_name = f"{gene_id}:ER{er_id}:EF{ef_id}"
        ea_data.append([gene_id, tx1_name, tx2_name, f"{tx1_name}|{tx2_name}", ef_name,
                        i.chrom, i.start, i.end, i.strand, 0, er_name, i.chrom, i.start,
                        i.end, i.strand])
        er_id += 1
    ea_df = pd.DataFrame(ea_data, columns=ea_df_cols)
    return ea_df

# !!! Below function potentially can replace format_identical_pair_ea (side = "none", or anything not "start", "end", or "both")
def format_fsm_pair_ea(tx1_bed_str, tx2_bed_str, tx1_name, tx2_name, gene_id, side_diff):
    ea_data = []
    er_id = 1
    tx1_bed_df = pd.DataFrame(tx1_bed_str,columns=['chrom','start','end','name','score','strand'])
    tx1_bed_df['start'] = tx1_bed_df['start'].astype(int)
    tx1_bed_df['end'] = tx1_bed_df['end'].astype(int)
    tx1_bed_df = tx1_bed_df.sort_values('start')
    tx2_bed_df = pd.DataFrame(tx2_bed_str,columns=['chrom','start','end','name','score','strand'])
    tx2_bed_df['start'] = tx2_bed_df['start'].astype(int)
    tx2_bed_df['end'] = tx2_bed_df['end'].astype(int)
    tx2_bed_df = tx2_bed_df.sort_values('start')
    max_er_id = len(tx1_bed_df)
    for index,i in tx1_bed_df.iterrows():
        ef_id = 1
        if er_id == 1 and (side_diff == "start" or side_diff == "both"):
            tx1_start = tx1_bed_df['start'].min()
            tx2_start = tx2_bed_df['start'].min()
            min_start = min(tx1_start,tx2_start)
            max_start = max(tx1_start,tx2_start)
            er_name = f"{gene_id}:ER{er_id}"
            ef_name = f"{gene_id}:ER{er_id}:EF{ef_id}"
            if min_start == tx1_start:
                tx_name = tx1_name
            else:
                tx_name = tx2_name
            ea_data.append([gene_id, tx1_name, tx2_name, f"{tx_name}", ef_name,
                            i.chrom, min_start, max_start, i.strand, 0, er_name, i.chrom, min_start,
                            i.end, i.strand])
            ef_id += 1
            ef_name = f"{gene_id}:ER{er_id}:EF{ef_id}"
            ea_data.append([gene_id, tx1_name, tx2_name, f"{tx1_name}|{tx2_name}", ef_name,
                            i.chrom, max_start, i.end, i.strand, 0, er_name, i.chrom, min_start,
                            i.end, i.strand])
        if er_id == max_er_id and (side_diff == "end" or side_diff == "both"):
            tx1_end = tx1_bed_df['end'].max()
            tx2_end = tx2_bed_df['end'].max()
            min_end = min(tx1_end,tx2_end)
            max_end = max(tx1_end,tx2_end)
            er_name = f"{gene_id}:ER{er_id}"
            ef_name = f"{gene_id}:ER{er_id}:EF{ef_id}"
            ea_data.append([gene_id, tx1_name, tx2_name, f"{tx1_name}|{tx2_name}", ef_name,
                            i.chrom, i.start, min_end, i.strand, 0, er_name, i.chrom, i.start,
                            max_end, i.strand])
            ef_id += 1
            ef_name = f"{gene_id}:ER{er_id}:EF{ef_id}"
            if max_end == tx1_end:
                tx_name = tx1_name
            else:
                tx_name = tx2_name
            ea_data.append([gene_id, tx1_name, tx2_name, f"{tx_name}", ef_name,
                            i.chrom, min_end, max_end, i.strand, 0, er_name, i.chrom, i.start,
                            max_end, i.strand])
        if (er_id != 1 and er_id != max_er_id) or (side_diff not in ["start","end","both"]):
            er_name = f"{gene_id}:ER{er_id}"
            ef_name = f"{gene_id}:ER{er_id}:EF{ef_id}"
            ea_data.append([gene_id, tx1_name, tx2_name, f"{tx1_name}|{tx2_name}", ef_name,
                            i.chrom, i.start, i.end, i.strand, 0, er_name, i.chrom, i.start,
                            i.end, i.strand])
        er_id += 1
    ea_df = pd.DataFrame(ea_data, columns=ea_df_cols)
    return ea_df

def do_ea(tx_data):
    """Perform event analysis using pybedtools on bed string data"""
    try:
        bed_data = prep_bed_for_ea(tx_data)
    except ValueError:
        raise
    if len(bed_data) == 4:
        ea_results, jct_catalog, transcript_distance = do_ea_pair(bed_data)
        return ea_results, jct_catalog, transcript_distance
    # full-gene, more transcripts than two
    else:
        ea_results, jct_catalog = do_ea_gene(bed_data)
        return ea_results, jct_catalog


def ea_pairwise(data, out_fhs, gene_id):
    "Do EA (Event Analysis) for a pair transcripts."
    transcripts = data.groupby("transcript_id")
    transcript_groups = transcripts.groups
    tx_data = {}
    for transcript in transcript_groups:
        transcript_df = data[data['transcript_id'] == transcript]
        tx_data[transcript] = transcript_df
    transcript_pairs = list(itertools.combinations(tx_data.keys(), 2))
    ea_df = pd.DataFrame(columns=ea_df_cols)
    jct_df = pd.DataFrame(columns=jct_df_cols)
    td_df = pd.DataFrame(columns=TD.td_df_cols)
    for tx_pair in transcript_pairs:
        # try:
        tx_df_1 = tx_data[tx_pair[0]]
        tx_df_2 = tx_data[tx_pair[1]]
        tx_pair_data = pd.concat([tx_df_1, tx_df_2])
        # try:
        ea_data, jct_data, td_data = do_ea(tx_pair_data)
        ea_df = ea_df.append(ea_data)
        jct_df = jct_df.append(jct_data)
        td_df = td_df.append(td_data)
        # except ValueError as e:
        #     logger.error(e)
        #     exit(1)
        #     #     # DEBUG, exit for now, just write to rejects later
        #     #     # continue
    return ea_df, jct_df, td_df


def process_single_file(infile, ea_mode, outdir, outfiles):
    """Compare all transcript pairs in a single GTF file."""
    logger.info("Input file: {}", infile)
    logger.debug("Output files: {}", outfiles)
    out_fhs = open_output_files(outdir, outfiles)
    out_fhs['ea_fh'].write_text(",".join(ea_df_cols) + '\n')
    out_fhs['jc_fh'].write_text(",".join(jct_df_cols) + '\n')
    data = read_exon_data_from_file(infile)
    genes = data.groupby("gene_id")
    transcripts = data.groupby(["gene_id", "transcript_id"])
    logger.info("Found {} genes and {} transcripts", len(genes), len(transcripts))

    for gene in genes.groups:
        gene_df = data[data['gene_id'] == gene]
        transcripts = gene_df.groupby("transcript_id")
        transcript_groups = transcripts.groups
        number_of_transcripts = len(transcript_groups)
        if number_of_transcripts == 1:
            logger.info("Gene {} has a single transcript. Skipping", gene)
            continue
        if ea_mode == 'gene':
            try:
                ea_data, jct_data = do_ea(gene_df)
                write_output(ea_data, out_fhs, 'ea_fh')
                write_output(jct_data, out_fhs, 'jc_fh')
            except ValueError as e:
                logger.error(e)
                continue
        else:
            out_fhs['td_fh'].write_text(",".join(TD.td_df_cols) + '\n')
            ea_data, jct_data, td_data = ea_pairwise(gene_df, out_fhs, gene)
            write_output(ea_data, out_fhs, 'ea_fh')
            write_output(jct_data, out_fhs, 'jc_fh')
            write_output(td_data, out_fhs, 'td_fh')


def add_suffix(value, suffix):
    """Add a column value suffix"""
    return value + suffix


def ea_two_files(f1_data, f2_data, out_fhs, gene_id):
    "Do EA (Event Analysis) for pairs of transcripts from two files for a gene."
    f1_transcripts = list(set(f1_data['transcript_id']))
    f2_transcripts = list(set(f2_data['transcript_id']))
    transcript_combos = list(itertools.product(f1_transcripts, f2_transcripts))
    logger.debug("Transcript combinations to process for {}: \n{}", gene_id, transcript_combos)
    ea_df = pd.DataFrame(columns=ea_df_cols)
    jct_df = pd.DataFrame(columns=jct_df_cols)
    td_df = pd.DataFrame(columns=TD.td_df_cols)
    for pair in transcript_combos:
        # try:
        tx_df_1 = f1_data[f1_data['transcript_id'] == pair[0]]
        tx_df_1_s = tx_df_1.assign(transcript_id=lambda x: x.transcript_id + '_d1')
        tx_df_2 = f2_data[f2_data['transcript_id'] == pair[1]]
        tx_df_2_s = tx_df_2.assign(transcript_id=lambda x: x.transcript_id + '_d2')
        tx_pair_data = pd.concat([tx_df_1_s, tx_df_2_s])
        try:
            ea_data, jct_data, td_data = do_ea(tx_pair_data)
            ea_df = ea_df.append(ea_data)
            jct_df = jct_df.append(jct_data)
            td_df = td_df.append(td_data)
        except ValueError:
            raise
        # except ValueError as e:
        #     logger.error(e)
        #     exit(1)
        #     #     # DEBUG, exit for now, just write to rejects later
        #     #     # continue
    return ea_df, jct_df


def process_two_files(infiles, outdir, outfiles):
    """Compare transcript pairs between two GTF files."""
    logger.info("Input files: {}", infiles)
    out_fhs = open_output_files(outdir, outfiles)
    out_fhs['ea_fh'].write_text(",".join(ea_df_cols) + '\n')
    out_fhs['jc_fh'].write_text(",".join(jct_df_cols) + '\n')
    out_fhs['td_fh'].write_text(",".join(TD.td_df_cols) + '\n')
    infile_1 = infiles[0]
    infile_2 = infiles[1]
    in_f1 = read_exon_data_from_file(infile_1)
    in_f2 = read_exon_data_from_file(infile_2)
    f1_gene_names = set(in_f1['gene_id'])
    f2_gene_names = set(in_f2['gene_id'])
    only_f1_genes = f1_gene_names.difference(f2_gene_names)
    only_f2_genes = f2_gene_names.difference(f1_gene_names)
    odd_genes = only_f1_genes.union(only_f2_genes)
    f1_odds = in_f1[in_f1['gene_id'].isin(only_f1_genes)]
    f2_odds = in_f2[in_f2['gene_id'].isin(only_f2_genes)]
    write_output(f1_odds, out_fhs, 'gtf1_fh')
    write_output(f2_odds, out_fhs, 'gtf2_fh')
    common_genes = f1_gene_names.difference(odd_genes)
    valid_f1 = in_f1[in_f1['gene_id'].isin(common_genes)]
    valid_f2 = in_f2[in_f2['gene_id'].isin(common_genes)]
    f1_genes = valid_f1.groupby("gene_id")
    f1_transcripts = valid_f1.groupby(["gene_id", "transcript_id"])
    f2_genes = valid_f2.groupby("gene_id")
    f2_transcripts = valid_f2.groupby(["gene_id", "transcript_id"])
    logger.info("Found {} genes and {} transcripts in {} file", len(f1_genes), len(f1_transcripts),
                infile_1)
    logger.info("Found {} genes and {} transcripts in {} file", len(f2_genes), len(f2_transcripts),
                infile_2)
    gene_list = list(set(valid_f1['gene_id']))
    logger.debug("Genes to process: \n{}", gene_list)
    for gene in gene_list:
        f1_data = valid_f1[valid_f1['gene_id'] == gene]
        f2_data = valid_f2[valid_f2['gene_id'] == gene]
        logger.debug(f1_data.head())
        logger.debug(f2_data.head())
        try:
            ea_data, jct_data, td_data = ea_two_files(f1_data, f2_data, out_fhs, gene)
        except ValueError as e:
            logger.error(e)
            continue
        write_output(ea_data, out_fhs, 'ea_fh')
        write_output(jct_data, out_fhs, 'jc_fh')
        write_output(td_data, out_fhs, 'td_fh')


def main():
    """Main function"""
    args = parse_args()
    setup_logging(args.debug, args.verbose, args.log_file)
    logger.debug("Args: {}", args)
    infiles = args.infiles
    outdir = handle_outdir(args)
    ea_mode = args.ea_mode
    if len(infiles) == 1:
        logger.debug("Single file {} analysis", ea_mode)
        outfiles = common_outfiles
        if ea_mode == "pairwise":
            outfiles.update(pairwise_outfiles)
        process_single_file(infiles[0], ea_mode, outdir, outfiles)
    else:
        logger.debug("Two files pairwise analysis")
        outfiles = common_outfiles
        outfiles.update(two_gtfs_outfiles)
        outfiles.update(pairwise_outfiles)
        process_two_files(infiles, outdir, outfiles)
    # The End


if __name__ == '__main__':
    main()
