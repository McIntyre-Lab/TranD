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
# from interval import interval
# from intervaltree import Interval
# from intervaltree import IntervalTree
from loguru import logger
from numpy import nan
from pathlib import Path
from pybedtools import BedTool


common_outfiles = {'ea_er_fh': 'gene_ea_exonic_regions.csv', 'ea_ef_fh':
                   'gene_ea_exonic_fragments.csv', 'ea_pairwise_fh': 'pairwise_ea.csv', 'td_fh':
                   'transcript_distances.csv'}

two_gtfs_outfiles = {'gtf1_fh': 'gtf1_only.gtf', 'gtf2_fh': 'gtf2_only.gtf', 'gtf_names_fh':
                     'gtf_metadata.csv'}


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
            logger.warning("EA mode is ignored for two GTF files - only pairwise EA is done.")
        else:
            logger.warning("No TD will be done as full-gene EA was specified.")
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


def setup_logging(debug, logfile):
    """Set the correct logging level and sinks."""
    logger.remove()
    level = logging.INFO
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


def open_output_files(outdir, outfiles):
    """Open all necessary output files and return filehandles"""
    out_fhs = {}
    try:
        for outfile in outfiles:
            logger.debug("Open fh: {}", outfiles[outfile])
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
    logger.debug("Drop unneeded columns: {}", drop_cols)
    logger.debug("Raw data top:\n{}", data.head())
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
    logger.debug("Exon data top:\n{}", data.head())
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


def write_ea_data(data, out_fhs):
    """Write results of event analysis to output files."""
    logger.debug("Writing EA data")
    data.to_csv(out_fhs['ea_pairwise_fh'], mode='a', header=False, index=False)


def prep_bed_for_ea(data):
    """Event analysis"""
    exons = {}
    logger.debug("Tx data for EA: {}", data)
    gene_id = data['gene_id'].unique()[0]
    logger.debug("GeneID: {}", gene_id)
    tx_names = list(data['transcript_id'].unique())
    logger.debug("Do EA on transcripts: {}", tx_names)
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
        logger.debug("Transcript data for {}:\n{}", tx, tx_data)
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


def do_ea_pair(data):
    """
    Event Analysis on a pair of transcripts
    """
    ea_df_cols = ['gene_id', 'transcript_1', 'transcript_2', 'transcript_id', 'ef_id', 'ef_chr',
                  'ef_start', 'ef_end', 'ef_strand', 'ef_ir_flag', 'er_id', 'er_chr', 'er_start',
                  'er_end', 'er_strand']
    ea_data = []
    gene_id = data['gene_id']
    tx_names = data['transcript_list']
    tx1_name, tx2_name = tx_names[0], tx_names[1]
    tx1_bed_str = data[tx1_name]
    tx1_bed = BedTool(tx1_bed_str).saveas()
    tx2_bed_str = data[tx2_name]
    tx2_bed = BedTool(tx2_bed_str).saveas()
    t2_t1_v_intersect = tx2_bed.intersect(tx1_bed, v=True)
    t1_t2_v_intersect = tx1_bed.intersect(tx2_bed, v=True)
    # Check if identical - enumerate EFs from tx1 if true and return
    if t1_t2_v_intersect == t2_t1_v_intersect:
        logger.info("ERs are identical for {}.", " / ".join(tx_names))
        # Represent the features as if EFs were generated from non-identical ERs to generalize
        er_id, ef_id = 1, 1
        logger.debug("ef_id,chrom,start,stop,strand")
        for i in tx1_bed:
            er_name = f"{gene_id}:ER{er_id}"
            ef_name = f"{gene_id}:ER{er_id}:EF{ef_id}"
            ea_data.append([gene_id, tx1_name, tx2_name, f"{tx1_name}|{tx2_name}", ef_name, i.chrom,
                            i.start, i.stop, i.strand, 0, er_name, i.chrom, i.start, i.stop,
                            i.strand])
        ea_df = pd.DataFrame(ea_data, columns=ea_df_cols)
        return ea_df
    # Transcripts are not identical
    else:
        strand = list(set([x.strand for x in tx1_bed]))[0]
        raw_ers = tx1_bed.cat(tx2_bed, postmerge=True).saveas()
        raw_ers_list = [str(x).split() for x in raw_ers]
        er_id = 1
        ers_list = []
        for i in raw_ers_list:
            i.extend([f"{gene_id}:ER{er_id}", '0', strand])
            ers_list.append("\t".join(i))
            er_id += 1
        ers_str = "\n".join(ers_list)
        ers_bed = BedTool(ers_str, from_string=True).saveas()
        logger.debug("ERs: \n{}", ers_bed)
        # Skipped exons - ERs and EFs are identical
        skipped_exons = []
        if len(t1_t2_v_intersect) > 0:
            for i in t1_t2_v_intersect:
                s_exon = i.name
                skipped_exons.append(s_exon)
                s_exon_bed = t1_t2_v_intersect.filter(lambda x: x.name == s_exon).saveas()
                logger.debug("Skipped exon: \n{}", s_exon_bed)
                exit("DEBUG")
        # Refactor from here onward since I found ERs apriori above in a new approach.
        # Don't write skipped exon data out, yet.
        # Keep all EFs in a data structure, enumerate sequentially after all EFs are described.
                ea_data.append([gene_id, tx1_name, tx2_name, tx1_name, '', i.chrom, i.start, i.stop,
                                i.strand, 0, '', i.chrom, i.start, i.stop, i.strand])
        if len(t2_t1_v_intersect) > 0:
            for i in t2_t1_v_intersect:
                skipped_exons.append(i.name)
                ea_data.append([gene_id, tx1_name, tx2_name, tx2_name, '', i.chrom, i.start, i.stop,
                                i.strand, 0, '', i.chrom, i.start, i.stop, i.strand])
        exit("DEBUG")
        # Remaining exons are overlapping, slice into individual EFs
        # To process overlapping exons filter out skipped exons
        # Filter produces an iterator. Do not consume prematurely or save to tmp file (performance?)
        tx1_o_bed = tx1_bed.filter(lambda x: x.name not in skipped_exons).saveas()
        logger.debug('TX1O: \n{}', tx1_o_bed)
        tx2_o_bed = tx2_bed.filter(lambda x: x.name not in skipped_exons).saveas()
        # Remove saveas after debugging
        raw_ers = tx1_o_bed.cat(tx2_o_bed, postmerge=True).saveas()
        logger.debug("ERs: \n{}", raw_ers)
        exit("DEBUG")

    out_df = pd.DataFrame(ea_data, columns=ea_df_cols)
    return out_df


def do_ea(tx_data):
    """Perform even analysis using pybedtools on bed string data"""
    bed_data = prep_bed_for_ea(tx_data)
    logger.debug("BED: {}", bed_data)
    if len(bed_data) == 4:
        ea_results = do_ea_pair(bed_data)
        return ea_results
    # full-gene, more transcripts than two
    else:
        raise NotImplementedError


def ea_pairwise(data, out_fhs, gene_id):
    "Do EA (Event Analysis) for a pair transcripts."
    logger.debug(f"EA for '{gene_id}' gene.")
    transcripts = data.groupby("transcript_id")
    transcript_groups = transcripts.groups
    tx_data = {}
    for transcript in transcript_groups:
        transcript_df = data[data['transcript_id'] == transcript]
        tx_data[transcript] = transcript_df
    transcript_pairs = list(itertools.combinations(tx_data.keys(), 2))
    logger.debug("All transcript pairs for {}: {}", gene_id, transcript_pairs)
    for tx_pair in transcript_pairs:
        try:
            tx_df_1 = tx_data[tx_pair[0]]
            tx_df_2 = tx_data[tx_pair[1]]
            tx_pair_data = pd.concat([tx_df_1, tx_df_2])
            try:
                ea_data = do_ea(tx_pair_data)
                return ea_data
            except ValueError as e:
                logger.error(e)
                exit(1)
                # DEBUG, exit now, write to rejects later
                # continue
        except ValueError as e:
            logger.error(e)


def process_single_file(infile, ea_mode, outdir, outfiles):
    """Compare all transcript pairs in a single GTF file."""
    logger.info("Input file: {}", infile)
    logger.debug("Output files: {}", outfiles)
    out_fhs = open_output_files(outdir, outfiles)
    data = read_exon_data_from_file(infile)
    # logger.debug("Data for EA:\n{}", data.head())
    genes = data.groupby("gene_id")
    transcripts = data.groupby(["gene_id", "transcript_id"])
    logger.info("Found {} genes and {} transcripts", len(genes), len(transcripts))
    for gene in genes.groups:
        gene_df = data[data['gene_id'] == gene]
        transcripts = gene_df.groupby("transcript_id")
        transcript_groups = transcripts.groups
        number_of_transcripts = len(transcript_groups)
        if number_of_transcripts == 1:
            logger.warning("Gene {} has a single transcript. Skipping", gene)
            continue
        if ea_mode == 'gene':
            try:
                ea_data = do_ea(gene_df)
                write_ea_data(ea_data, out_fhs)
            except ValueError as e:
                logger.error(e)
                continue
        else:
            ea_data = ea_pairwise(gene_df, out_fhs, gene)
            write_ea_data(ea_data, out_fhs)


def process_two_files(infiles):
    """Compare transcript pairs between two GTF files."""
    logger.info("Input files: {}", infiles)
    exit(logger.error("Two file comparison is not implemented, yet"))


def main():
    """Main function"""
    args = parse_args()
    setup_logging(args.debug, args.log_file)
    logger.debug("Args: {}", args)
    infiles = args.infiles
    outdir = handle_outdir(args)
    ea_mode = args.ea_mode
    if len(infiles) == 1:
        outfiles = common_outfiles
        process_single_file(infiles[0], ea_mode, outdir, outfiles)
    else:
        outfiles = common_outfiles
        outfiles.update(two_gtfs_outfiles)
        process_two_files(infiles, outdir, outfiles)
    # The End


if __name__ == '__main__':
    main()
