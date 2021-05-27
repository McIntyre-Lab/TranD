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
from dataclasses import dataclass
# from itertools import chain
from loguru import logger
from numpy import nan
from pathlib import Path
from pybedtools import BedTool
# from pybedtools import Interval


# CONFIGURATION
common_outfiles = {'ea_fh': 'event_analysis.csv', 'jc_fh': 'junction_catalog.csv', 'er_fh':
                   'event_analysis_er.csv', 'ef_fh': 'event_analysis_ef.csv'}

ea_df_cols = ['gene_id', 'transcript_1', 'transcript_2', 'transcript_id', 'ef_id', 'ef_chr',
              'ef_start', 'ef_end', 'ef_strand', 'ef_ir_flag', 'er_id', 'er_chr', 'er_start',
              'er_end', 'er_strand']
jct_df_cols = ['gene_id', 'transcript_id', 'coords']

er_df_cols = ['gene_id', 'er_id', 'er_chr', 'er_start', 'er_end', 'er_strand', 'er_exon_ids',
              'er_transcript_ids', 'gene_transcript_ids', 'exons_per_er', 'transcripts_per_er',
              'transcripts_per_gene',  'er_ir_flag', 'er_annotation_frequency']

ef_df_cols = ['gene_id', 'er_id', 'ef_id', 'ef_chr', 'ef_start', 'ef_end', 'ef_strand',
              'ef_exon_ids', 'ef_transcript_ids', 'exons_per_ef', 'transcripts_per_ef',
              'ef_ir_flag', 'ea_annotation_frequency']


# Later when TD has been added
# common_outfiles = {'ea_er_fh': 'gene_ea_exonic_regions.csv', 'ea_ef_fh':
#                    'gene_ea_exonic_fragments.csv', 'ea_fh': 'pairwise_ea.csv', 'td_fh':
#                    'transcript_distances.csv'}

two_gtfs_outfiles = {'gtf1_fh': 'gtf1_only.gtf', 'gtf2_fh': 'gtf2_only.gtf'}
# , 'gtf_names_fh': 'gtf_metadata.csv'}

KEEP_IR = False


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
        "-k", "--keepir",
        action="store_true",
        help="Keep transcripts with Intron Retention events in full-gene EA. Default: remove",
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


# Data structures for ERs and EFs. Use a mutable dataclass, so we could add transcript names and
# counts during EA
@dataclass
class ER:
    gene_id: str
    id: str
    chrom: str
    start: int
    end: int
    strand: str
    ex_ids: list
    tx_ids: list
    gene_tx_ids: list
    ex_num: int = 0
    tx_num: int = 0
    gene_tx_num: int = 0
    ir_flag: int = 0
    er_freq: str = ''


@dataclass
class EF:
    gene_id: str
    er_id: str
    ef_id: str
    chrom: str
    start: int
    end: int
    strand: str
    ex_ids: list
    tx_ids: list
    ex_num: int = 0
    tx_num: int = 0
    ir_flag: int = 0
    ef_freq: str = ''


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


def remove_ir_transcripts(tx_data, introns):
    """Remove transcripts that contain Intron Retention events i.e. an
    intron is encompassed by an exon
    """
    to_remove = []
    for e_tx in tx_data:
        REMOVED = False
        e_ints = [(int(x.start), int(x.end)) for x in tx_data[e_tx]]
        e_ints.sort(key=lambda x: x[0])
        for i in introns:
            if REMOVED:
                break
            for e in e_ints:
                # This is the exon this intron may overlap with
                if i[0] >= e[0] and i[0] <= e[1]:
                    # Intron is fully encompassed within the exon
                    if i[1] <= e[1]:
                        logger.debug("IR detected, removing {}", e_tx)
                        # logger.debug("Exon: {}:{}, intron: {}:{}", e[0], e[1], i[0], i[1])
                        to_remove.append(e_tx)
                        REMOVED = True
                        break
                    else:
                        # Stop processing this transcript's exons, move on to next intron
                        continue
    for ir_tx in to_remove:
        del tx_data[ir_tx]
    return tx_data


def do_ea_gene(data):
    """
    Event Analysis on a pair of transcripts
    """
    gene_id = data['gene_id']
    logger.debug("Gene EA on {}", gene_id)
    # logger.debug(data)
    tx_names = data['transcript_list']
    logger.debug("Transcripts in {}: {}", gene_id, tx_names)
    tx_data = {}
    tx_coords = {}
    all_introns = []
    for tx_name in tx_names:
        tx_bed_str = data[tx_name]
        starts = [int(x[1]) for x in tx_bed_str]
        ends = [int(x[2]) for x in tx_bed_str]
        chrom = data[tx_name][0][0]
        strand = data[tx_name][0][5]
        feature = f"{tx_name}_transcript"
        left_edge = min(starts)
        right_edge = max(ends)
        tx_exons = BedTool(tx_bed_str).saveas()
        tx_data[tx_name] = tx_exons
        tx_coords[tx_name] = [left_edge, right_edge]
        if not KEEP_IR:
            tx_full = BedTool(f"{chrom} {left_edge} {right_edge} {feature} 0 {strand}",
                              from_string=True)
            tx_introns = tx_full.subtract(tx_exons)
            i_coords = [(int(x.start), int(x.end)) for x in tx_introns]
            all_introns.extend(i_coords)
        # logger.debug("Gene EA Raw Tx: {}\n{}", tx_name, tx_bed_str)
    intron_data = list(set(all_introns))
    intron_data.sort(key=lambda x: x[0])
    # Removal of IR containing transcripts
    if not KEEP_IR:
        tx_data = remove_ir_transcripts(tx_data, intron_data)
        old_tx_names = tx_names
        tx_names = tx_data.keys()
        tx_coords_delete = [k for k in old_tx_names if k not in tx_names]
        for k in tx_coords_delete:
            del tx_coords[k]
    else:
        logger.info("Not attempting to remove IR-containing transcripts")
    # Junction Catalog creation
    junction_data = []
    for tx_name in tx_names:
        tx_jc_data = create_junction_catalog(gene_id, tx_name, tx_data[tx_name])
        junction_data.extend(tx_jc_data)
        for jc in tx_jc_data:
            tx_coords[tx_name].append(jc[2])
    junction_df = pd.DataFrame(junction_data, columns=jct_df_cols)
    # Event Analysis
    er_data, ef_data = ea_analysis(gene_id, tx_data, tx_coords)
    er_df = pd.DataFrame(er_data, columns=er_df_cols)
    ef_df = pd.DataFrame(ef_data, columns=ef_df_cols)
    return er_df, ef_df, junction_df


def do_ea_pair(data):
    """
    Event Analysis on a pair of transcripts
    """
    gene_id = data['gene_id']
    tx_names = data['transcript_list']
    tx1_name, tx2_name = tx_names[0], tx_names[1]
    tx1_bed_str = data[tx1_name]
    tx1_bed = BedTool(tx1_bed_str).saveas()
    tx2_bed_str = data[tx2_name]
    tx2_bed = BedTool(tx2_bed_str).saveas()
    logger.debug("Comparing {} vs {}", tx1_name, tx2_name)
    junction_data = create_junction_catalog(gene_id, tx1_name, tx1_bed)
    junction_data.extend(create_junction_catalog(gene_id, tx2_name, tx2_bed))
    # Check if identical - enumerate EFs from tx1 if true and return
    t2_t1_v_intersect = tx2_bed.intersect(tx1_bed, v=True)
    t1_t2_v_intersect = tx1_bed.intersect(tx2_bed, v=True)
    # Double step verification to get a bit more performance
    if t1_t2_v_intersect == t2_t1_v_intersect:
        sub_2_from_1 = tx2_bed.subtract(tx1_bed)
        sub_1_from_2 = tx1_bed.subtract(tx2_bed)
        if sub_2_from_1 == sub_1_from_2:
            # Transcripts are identical
            logger.info("ERs are identical for {}.", " / ".join(tx_names))
            ea_data = format_identical_pair_ea(tx1_bed, tx2_bed, tx1_name, tx2_name, gene_id)
        else:
            # Transcripts are not identical
            ea_data = er_ea_analysis(tx1_bed, tx2_bed, tx1_name, tx2_name, gene_id)
    # Transcripts are not identical
    else:
        ea_data = er_ea_analysis(tx1_bed, tx2_bed, tx1_name, tx2_name, gene_id)
    out_df = pd.DataFrame(ea_data, columns=ea_df_cols)
    junction_df = pd.DataFrame(junction_data, columns=jct_df_cols)
    return out_df, junction_df


def er_ea_analysis(tx1_bed, tx2_bed, tx1_name, tx2_name, gene_id):
    """
    Generate ERs (Exonic Regions) and EFs (Exonic Fragments) and analyze events in a pair of
    transcripts.
    """
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


def ea_analysis(gene_id, tx_data, tx_coords):
    """
    Generate ERs (Exonic Regions) and EFs (Exonic Fragments) and analyze events in transcripts.
    Generalize to do both full-gene and transcript pairs to replace er_ea_analysis.

    Required outputs:
    gene_id
    gene_transcript_ids
    transcripts_per_gene
    er_id
    er_chr
    er_start
    er_end
    er_strand
    er_flag_ir ( zero if IRs were removed)
    exons_per_er
    er_exon_ids = Piped ("|") list of exon IDs that are contained within the exon region
    transcripts_per_er
    er_transcript_ids = Piped ("|") list of transcript IDs the exon region is present in
    er_annotation_frequency = “unique”, “common”, “constitutive” (one, many, all) per txs
    ef_id
    ef_chr
    ef_start
    ef_end
    ef_strand
    ef_flag_ir (zero if IRs were removed)
    ef_exon_ids
    exons_per_ef
    transcripts_per_ef
    ef_transcript_ids = Piped ("|") list of transcript IDs the exon fragment is present in
    ef_annotation_frequency = “unique”, “common”, “constitutive” (one, many, all) per txs
    """
    # Use the junction catalog and transcripts start/end to check for identical transcripts
    logger.debug("Performing EA analysis for {} gene", gene_id)
    logger.debug("TX Names: {}", list(tx_data.keys()))
    # logger.debug("TX Coords: {}", tx_coords)
    gene_transcript_ids = list(tx_data.keys())
    total_number_of_transcripts = len(gene_transcript_ids)
    tx_ident_check = {}
    for i in tx_coords:
        tx_structure_str_list = [str(x) for x in tx_coords[i]]
        tx_structure = ",".join(tx_structure_str_list)
        if tx_structure in tx_ident_check:
            tx_ident_check[tx_structure].append(i)
        else:
            tx_ident_check[tx_structure] = [i]
    identical_tx_list = [x for x in tx_ident_check.values() if len(x) > 1]
    identical_transcripts = {x[0]: x[1:] for x in identical_tx_list}
    logger.debug("Identical records: {}", identical_transcripts)
    # Stash names for EA output, but remove duplicate data
    duplicate_txs = [x[1:] for x in identical_tx_list]
    logger.debug("Duplicate records: {}", duplicate_txs)
    for dupe_list in duplicate_txs:
        for item in dupe_list:
            del tx_data[item]
    logger.debug("Unique TX Data: {}", list(tx_data.keys()))
    # Combine transcripts and generate ERs
    strands = []
    all_tx = None
    for i in tx_data:
        tx = tx_data[i]
        if not all_tx:
            all_tx = tx
        else:
            all_tx = all_tx.cat(tx, postmerge=True)
        strand = list(set([x.strand for x in tx]))[0]
        strands.append(strand)
    strand = list(set(strands))[0]
    # Is this a valid check?
    if len(strand) > 1:
        raise ValueError("Multiple strands detected when there should be one")
    logger.debug("Strand: '{}'", strand)
    raw_ers_list = [str(x).split() for x in all_tx]
    er_id = 1
    ers_list = []
    for i in raw_ers_list:
        i.extend([f"{gene_id}:ER{er_id}", '0', strand])
        ers_list.append("\t".join(i))
        er_id += 1
    ers_str = "\n".join(ers_list)
    ers_bed = BedTool(ers_str, from_string=True)
    er_data = {}
    # Coordinates of ER ends to check when creating EFs
    er_ends = []
    er_range_sets = {}
    for er in ers_bed:
        er_data[er.name] = ER(gene_id=gene_id, id=er.name, chrom=er.chrom, start=er.start,
                              end=er.end, strand=er.strand, ex_ids=[], tx_ids=[],
                              gene_tx_ids=gene_transcript_ids)
        er_ends.append(er.end)
        er_range_sets[er.name] = set(range(er.start, er.end + 1))
    # Generate EFs
    # Catalog starting and ending points (0-coordinate as we go from left side only)
    sts_ends = []
    tx_ranges = {}
    for i in tx_data:
        tx = tx_data[i]
        for ex in tx:
            # logger.debug(f"{ex.name}: {ex.start}/{ex.end}")
            # There are many ways to check if exon start is in ER, pick one at random
            if i in tx_ranges:
                tx_ranges[ex.name].append(ex.start, ex.end, i)
            else:
                tx_ranges[ex.name] = (ex.start, ex.end, i)
            if ex.start not in sts_ends:
                sts_ends.append(ex.start)
            if (ex.end) not in sts_ends:
                sts_ends.append(ex.end)
    sts_ends.sort()
    tx_range_sets = {}
    tx_exon_map = {}
    for i in tx_ranges:
        tx = tx_ranges[i]
        tx_range_set = set(range(tx[0], tx[1] + 1))
        tx_range_sets[i] = tx_range_set
        tx_exon_map[i] = tx[2]
    # Define EFs
    ef_list = []
    ef_start = sts_ends[0]
    ef_end = None
    for i in sts_ends[1:]:
        # Start of an ER
        if not ef_start and not ef_end:
            ef_start = i
            continue
        ef_end = i
        ef_list.append((ef_start, ef_end))
        if ef_end in er_ends:
            # Skip the intron
            ef_start = None
            ef_end = None
            continue
        else:
            ef_start = ef_end
    # Process all ERs and EFs vs transcript exons/transcripts to generate requested EA output
    # Sets are in 'er_range_sets' and 'tx_range_sets'
    for er_id in er_data:
        for ex_id in tx_exon_map:
            exon_data = tx_ranges[ex_id]
            # Check if exon start is in the ER
            if exon_data[0] in er_range_sets[er_id]:
                # Exon is in the ER, record exon, tx, and duplicates if any
                tx_id = exon_data[2]
                if ex_id not in er_data[er_id].ex_ids:
                    er_data[er_id].ex_ids.append(ex_id)
                if tx_id not in er_data[er_id].tx_ids:
                    er_data[er_id].tx_ids.append(tx_id)
                # Add duplicate txs and their exons if exist
                if tx_id in identical_transcripts:
                    duplicate_tx_list = identical_transcripts[tx_id]
                    for duplicate_tx in duplicate_tx_list:
                        if duplicate_tx not in er_data[er_id].tx_ids:
                            er_data[er_id].tx_ids.append(duplicate_tx)
                        ex_number = ex_id.split('_exon_')[1]
                        duplicate_ex = f"{duplicate_tx}_exon_{ex_number}"
                        if duplicate_ex not in er_data[er_id].ex_ids:
                            er_data[er_id].ex_ids.append(ex_id)
    for er in er_data:
        ex_num = len(er_data[er].ex_ids)
        tx_num = len(er_data[er].tx_ids)
        er_data[er].ex_num = ex_num
        er_data[er].tx_num = tx_num
        if tx_num == 1:
            er_data[er].er_freq = 'unique'
        elif tx_num == total_number_of_transcripts:
            er_data[er].er_freq = 'constitutive'
        else:
            er_data[er].er_freq = 'common'

    # Process EFs vs Exons and Transcripts
    ef_data = {}
    ef_start_set = set([x[0] for x in ef_list])
    ef_start_end_map = {x[0]: x[1] for x in ef_list}
    er_ef_map = {}
    # Match agains ERs to get id and common data
    for er in er_range_sets:
        er_set = er_range_sets[er]
        er_ef_starts = list(er_set.intersection(ef_start_set))
        for ef_start in er_ef_starts:
            er_ef_map[er] = er_ef_starts
    for er_id in er_ef_map:
        ef_id_ord = 1
        for ef_start in er_ef_map[er_id]:
            ef_id = f"{er_id}:EF{ef_id_ord}"
            ef_id_ord += 1
            ef_data[ef_id] = EF(gene_id=gene_id, er_id=er, ef_id=ef_id, chrom=er_data[er_id].chrom,
                                start=ef_start, end=ef_start_end_map[ef_start],
                                strand=er_data[er_id].strand, ex_ids=[], tx_ids=[])
    # Match against exons
    for ex_id in tx_ranges:
        ex_set = set(range(tx_ranges[ex_id][0], tx_ranges[ex_id][1] + 1))
        for ef_id in ef_data:
            ef_start = ef_data[ef_id].start
            if ef_start in ex_set:
                # EF is a parf of an exon - are there any duplicate exons and transcripts?
                tx_id = tx_exon_map[ex_id]
                if ex_id not in ef_data[ef_id].ex_ids:
                    ef_data[ef_id].ex_ids.append(ex_id)
                if tx_id not in ef_data[ef_id].tx_ids:
                    ef_data[ef_id].tx_ids.append(tx_id)
                # Add duplicate txs and their exons if exist
                if tx_id in identical_transcripts:
                    duplicate_tx_ids = identical_transcripts[tx_id]
                    ex_number = ex_id.split('_exon_')[1]
                    for dupe_tx_id in duplicate_tx_ids:
                        if dupe_tx_id not in ef_data[ef_id].tx_ids:
                            ef_data[ef_id].tx_ids.append(dupe_tx_id)
                        dupe_ex_id = f"{dupe_tx_id}_exon_{ex_number}"
                        if dupe_ex_id not in ef_data[ef_id].ex_ids:
                            ef_data[ef_id].ex_ids.append(dupe_ex_id)
    for ef in ef_data:
        ex_num = len(ef_data[ef].ex_ids)
        tx_num = len(ef_data[ef].tx_ids)
        ef_data[ef].ex_num = ex_num
        ef_data[ef].tx_num = tx_num
        if tx_num == 1:
            ef_data[ef].ef_freq = 'unique'
        elif tx_num == total_number_of_transcripts:
            ef_data[ef].ef_freq = 'constitutive'
        else:
            ef_data[ef].ef_freq = 'common'
    # er_df_cols = ['gene_id', 'er_id', 'er_chr', 'er_start', 'er_end', 'er_strand', 'er_exon_ids',
    #         'er_transcript_ids', 'gene_transcript_ids', 'exons_per_er', 'transcripts_per_er',
    #         'transcripts_per_gene',  'er_ir_flag', 'er_annotation_frequency']
    er_out_data = []
    for er in er_data:
        i = er_data[er]
        er_out_data.append([i.gene_id, i.id, i.chrom, i.start, i.end, i.strand, "|".join(i.ex_ids),
                            "|".join(i.tx_ids), "|".join(i.gene_tx_ids), i.ex_num, i.tx_num,
                            i.gene_tx_num, i.ir_flag, i.er_freq])
    # ef_df_cols = ['gene_id', 'er_id', 'ef_id', 'ef_chr', 'ef_start', 'ef_end', 'ef_strand',
    #         'ef_exon_ids', 'ef_transcript_ids', 'exons_per_ef', 'transcripts_per_ef',
    #         'ef_ir_flag', 'ea_annotation_frequence']
    ef_out_data = []
    for ef in ef_data:
        i = ef_data[ef]
        ef_out_data.append([i.gene_id, i.er_id, i.ef_id, i.chrom, i.start, i.end, i.strand,
                            "|".join(i.ex_ids), "|".join(i.tx_ids), i.ex_num, i.tx_num, i.ir_flag,
                            i.ef_freq])
    return er_out_data, ef_out_data


def do_ea(tx_data):
    """Perform even analysis using pybedtools on bed string data"""
    try:
        bed_data = prep_bed_for_ea(tx_data)
    except ValueError:
        raise
    if len(bed_data) == 4:
        ea_results, jct_catalog = do_ea_pair(bed_data)
        return ea_results, jct_catalog
    else:
        # full-gene, all transcripts
        er_df, ef_df, jct_catalog = do_ea_gene(bed_data)
        return er_df, ef_df, jct_catalog


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
    for tx_pair in transcript_pairs:
        # try:
        tx_df_1 = tx_data[tx_pair[0]]
        tx_df_2 = tx_data[tx_pair[1]]
        tx_pair_data = pd.concat([tx_df_1, tx_df_2])
        # try:
        ea_data, jct_data = do_ea(tx_pair_data)
        ea_df = ea_df.append(ea_data)
        jct_df = jct_df.append(jct_data)
        # except ValueError as e:
        #     logger.error(e)
        #     exit(1)
        #     #     # DEBUG, exit for now, just write to rejects later
        #     #     # continue
    return ea_df, jct_df


def process_single_file(infile, ea_mode, outdir, outfiles):
    """Compare all transcript pairs in a single GTF file."""
    logger.info("Input file: {}", infile)
    logger.debug("Output files: {}", outfiles)
    out_fhs = open_output_files(outdir, outfiles)
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
                out_fhs['ea_fh'].unlink(missing_ok=True)
                er_data, ef_data, jct_data = do_ea(gene_df)
                logger.debug("Writing data out to dis: {}", out_fhs)
                out_fhs['er_fh'].write_text(",".join(er_df_cols) + '\n')
                out_fhs['ef_fh'].write_text(",".join(ef_df_cols) + '\n')
                write_output(er_data, out_fhs, 'er_fh')
                write_output(ef_data, out_fhs, 'ef_fh')
                write_output(jct_data, out_fhs, 'jc_fh')
                exit("DEBUG")
            except ValueError as e:
                logger.error(e)
                continue
        else:
            out_fhs['ea_fh'].write_text(",".join(ea_df_cols) + '\n')
            ea_data, jct_data = ea_pairwise(gene_df, out_fhs, gene)
            write_output(ea_data, out_fhs, 'ea_fh')
            write_output(jct_data, out_fhs, 'jc_fh')


def ea_two_files(f1_data, f2_data, out_fhs, gene_id):
    "Do EA (Event Analysis) for pairs of transcripts from two files for a gene."
    f1_transcripts = list(set(f1_data['transcript_id']))
    f2_transcripts = list(set(f2_data['transcript_id']))
    transcript_combos = list(itertools.product(f1_transcripts, f2_transcripts))
    logger.debug("Transcript combinations to process for {}: \n{}", gene_id, transcript_combos)
    ea_df = pd.DataFrame(columns=ea_df_cols)
    jct_df = pd.DataFrame(columns=jct_df_cols)
    for pair in transcript_combos:
        # try:
        tx_df_1 = f1_data[f1_data['transcript_id'] == pair[0]]
        tx_df_1_s = tx_df_1.assign(transcript_id=lambda x: x.transcript_id + '_d1')
        tx_df_2 = f2_data[f2_data['transcript_id'] == pair[1]]
        tx_df_2_s = tx_df_2.assign(transcript_id=lambda x: x.transcript_id + '_d2')
        tx_pair_data = pd.concat([tx_df_1_s, tx_df_2_s])
        try:
            ea_data, jct_data = do_ea(tx_pair_data)
            ea_df = ea_df.append(ea_data)
            jct_df = jct_df.append(jct_data)
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
            ea_data, jct_data = ea_two_files(f1_data, f2_data, out_fhs, gene)
        except ValueError as e:
            logger.error(e)
            continue
        write_output(ea_data, out_fhs, 'ea_fh')
        write_output(jct_data, out_fhs, 'jc_fh')


def main():
    """Main function"""
    args = parse_args()
    setup_logging(args.debug, args.verbose, args.log_file)
    logger.debug("Args: {}", args)
    infiles = args.infiles
    outdir = handle_outdir(args)
    ea_mode = args.ea_mode
    if len(infiles) == 1:
        logger.debug("Single file pairwise analysis")
        outfiles = common_outfiles
        if args.keepir:
            global KEEP_IR
            KEEP_IR = True
        process_single_file(infiles[0], ea_mode, outdir, outfiles)
    else:
        logger.debug("Two files pairwise analysis")
        outfiles = common_outfiles
        outfiles.update(two_gtfs_outfiles)
        process_two_files(infiles, outdir, outfiles)
    # The End


if __name__ == '__main__':
    main()
