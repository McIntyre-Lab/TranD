#! /usr/bin/env python3
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2021 Oleksandr Moskalenko <om@rc.ufl.edu>
#
# Distributed under terms of the MIT license.

"""
Main module for the TranD package. Define a CLI and drive execution of all analyses.
"""

import argparse
import csv
import logging
import re
import os
import sys
from loguru import logger
from numpy import nan
import pandas as pd
from pathlib import Path
from pybedtools import cleanup
from trand.event_analysis import process_single_file
from trand.event_analysis import process_two_files


# CONFIGURATION
# Output file selections
# One file:
# common, gene OR pairwise
# Two files:
# common, two_gtfs, pairwise
common_outfiles = {
    "ea_fh": "event_analysis.csv",
    "jc_fh": "junction_catalog.csv",
    "er_fh": "event_analysis_er.csv",
    "ef_fh": "event_analysis_ef.csv",
}
pairwise_outfiles = {"td_fh": "pairwise_transcript_distance.csv"}
gene_outfiles = {"er_fh": "event_analysis_er.csv", "ef_fh": "event_analysis_ef.csv"}
two_gtfs_outfiles = {
    "gtf1_fh": "gtf1_only.gtf",
    "gtf2_fh": "gtf2_only.gtf",
    "md_fh": "minimum_pairwise_transcript_distance.csv",
}
consol_outfiles = {'key_fh': 'transcript_id_2_consolidation_id.csv', 'consol_gtf_fh':
                   'consolidated_transcriptome.gtf'}


def parse_args(print_help=False):
    """Parse command-line arguments"""

    class MyParser(argparse.ArgumentParser):
        """Subclass ArgumentParser for better help printing"""

        def error(self, message):
            sys.stderr.write("error: %s\n" % message)
            self.print_help()
            sys.exit(2)

    parser = MyParser(
        description="Perform transcript distance, complexity and "
                    "transcriptome comparison analyses."
    )
    parser.add_argument(
        dest="infiles",
        metavar="input_file",
        type=str,
        nargs="+",
        help="One or two input GTF file(s).",
    )
    parser.add_argument(
        "-o",
        "--outdir",
        action="store",
        type=str,
        required=False,
        help="Output directory, created if missing. Default: current directory.",
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
        "--consolidate",
        dest="consolidate",
        action="store_true",
        help="""Used with 1 GTF input file. Consolidate transcripts remove 5'/3' transcript end
        variation in redundantly spliced transcripts) with identical junctions prior to complexity
        calculations, events and summary plotting. Default: No consolidation""",
    )
    parser.add_argument(
        "--consolPrefix",
        dest="consol_prefix",
        type=str,
        default="tr",
        help="""Used with 1 GTF input file. Requires '--consolidate' flag. Specify the prefix to use
        for consolidated transcript_id values. Prefix must be alphanumeric with no spaces.
        Underscore (\"_\") is the only allowed special character. Default: 'tr'""",
    )
    parser.add_argument(
        "-c",
        "--complexityOnly",
        dest="complexity_only",
        action="store_true",
        help="""Used with 1 or 2 GTF input file(s). Output only complexity measures. If used in
        presence of the '--consolidate' flag, complexity is calculated on the consolidated GTF(s).
        Default: Perform all analyses and comparisons including complexity calculations""",
    )
    parser.add_argument(
        "-e",
        "--ea",
        dest="ea_mode",
        type=str,
        choices=["pairwise", "gene"],
        default="pairwise",
        help="""Specify type of within gene transcript comparison: pairwise - Used with 1 or 2 GTF
        input files. Compare pairs of transcripts within a gene. gene - Used iwth 1 GTF input file.
        Compare all transcripts within a gene Default: pairwise""",
    )
    parser.add_argument(
        "-k",
        "--keepir",
        dest="keep_ir",
        action="store_true",
        help="""Keep transcripts with Intron Retention(s) when generating transcript events. Only
        used with 1 GTF input file. Default: remove""",
    )
    parser.add_argument(
        "-p",
        "--pairs",
        type=str,
        choices=["all", "both", "first", "second"],
        dest="out_pairs",
        default="both",
        help="""Used with 2 GTF input files. The TranD metrics can be for all transcript pairs in
        both GTF files or for a subset of transcript pairs using the following options: both -
        Trand metrics for the minimum pairs in both GTF files, first - TranD metrics for the
        minimum pairs in the first GTF file, second - TranD metrics for the minimum pairs in the
        second GTF file all - TranD metrics for all transcript pairs in both GTF files Default:
        both""",
    )
    parser.add_argument(
        "-1",
        "--name1",
        dest="name1",
        default="d1",
        required=False,
        help="""Used with 2 GTF input files. User-specified name to be used for labeling output
        files related to the first GTF file. Name must be alphanumeric, can only include \"_\"
        special character and not contain any spaces. Default: d1""",
    )
    parser.add_argument(
        "-2",
        "--name2",
        dest="name2",
        default="d2",
        required=False,
        help="""Used with 2 GTF input files. User-specified name to be used for labeling output
        files related to the second GTF file. Name must be alphanumeric, can only include \"_\"
        special character and not contain any spaces. Default: d2""",
    )
    parser.add_argument(
        "-n",
        "--cpus",
        dest="cpu_cores",
        type=int,
        default=1,
        required=False,
        help="Number of CPU cores to use for parallelization. Default: 1",
    )

    parser.add_argument(
        "-f",
        "--force",
        action="store_true",
        help="Force overwrite existing output directory and files within.",
    )
    parser.add_argument(
        "-s",
        "--skip-plots",
        dest="skip_plots",
        action="store_true",
        help="Skip generation of all plots.",
    )
    parser.add_argument(
        "-i",
        "--skip-intermediate",
        dest="skip_interm",
        action="store_true",
        help="Skip intermediate file output (junction and exon region/fragment files).",
    )
    parser.add_argument("-v", "--verbose", action="store_true", help="Verbose output")
    parser.add_argument("-d", "--debug", action="store_true", help=argparse.SUPPRESS)
    if print_help:
        parser.print_help()
        sys.exit(0)
    args = parser.parse_args()
    if len(args.infiles) > 2:
        print("\nToo many input files - pass one or two GTF/GFF files as input.\n")
        parser.print_help()
        sys.exit(2)
    if args.ea_mode == "gene":
        if len(args.infiles) > 1:
            logger.warning(
                "EA 'gene' mode is ignored for two GTF files - only pairwise is done."
            )
    regex = re.compile(r'^\w+$', re.ASCII)
    # Validate prefixes
    if not regex.match(args.consol_prefix):
        logger.error("Invalid prefix format for consolidated transcript_id values: "
                     "Must be alphanumeric."
                     "Only '_' (underscore) special character is allowed"
                     )
        parser.print_help()
        sys.exit(2)
    if not regex.match(args.name1):
        logger.error(
            "Invalid name for dataset 1: Must be alphanumeric and can only "
            "include '_' special character"
        )
        parser.print_help()
        sys.exit(2)
    if not regex.match(args.name2):
        logger.error(
            "Invalid name for dataset 2: Must be alphanumeric and can only "
            "include '_' special character"
        )
        parser.print_help()
        sys.exit(2)
    # Multiprocessing checks
    if args.cpu_cores < 1:
        logger.error(
            "Invalid value for the number of CPU cores. Must be 1 or greater."
        )
        parser.print_help()
        sys.exit(2)
    # os.sched_getaffinity is accurate and linux HPC specific. os.cpu_count = total system cores.
    try:
        avail_cores = len(os.sched_getaffinity(0))
    except AttributeError:
        avail_cores = os.cpu_count()
    if args.cpu_cores > avail_cores:
        # Does this have to be an error since it leads to lowered performance?
        logger.warning("Requested CPU cores exceed the number of available cores!")
    return args


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


def get_gtf_attribute(transcript_id, gene_id):
    return f'transcript_id "{transcript_id}"; gene_id "{gene_id}";'


def prepare_outdir(args):
    if not args.outdir:
        outdir = Path.cwd()
    else:
        outdir = Path(args.outdir)
        logger.debug("Output directory: {}", str(outdir))
        if outdir.exists():
            if not args.force:
                exit("Not overwriting existing output directory without -f|--force. Exiting.")
        else:
            outdir.mkdir(parents=True, exist_ok=True)


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


def write_output(data, out_fhs, fh_name):
    """Write results of event analysis to output files."""
    data.to_csv(out_fhs[fh_name], mode='a', header=False, index=False)


def write_gtf(data, out_fhs, fh_name):
    """Write output gtf files."""
    if len(data) == 0:
        return
    else:
        data.loc[:, 'source'] = "TranD"
        data.loc[:, 'feature'] = "exon"
        data.loc[:, 'score'] = "."
        data.loc[:, 'frame'] = "."
        data.loc[:, 'attribute'] = data.apply(lambda x: get_gtf_attribute(x['transcript_id'],
                                              x['gene_id']), axis=1)
        output_column_names = ['seqname', 'source', 'feature', 'start', 'end', 'score', 'strand',
                               'frame', 'attribute']
        data = data.reindex(columns=output_column_names)
        data.to_csv(out_fhs[fh_name], sep="\t", mode='a', index=False, header=False,
                    doublequote=False, quoting=csv.QUOTE_NONE)


def read_exon_data_from_file(infile):
    """
    Create a pandas dataframe with exon records from a gtf file
    Raw gtf data:
    seqname source  feature  start end  score strand frame attributes
    2L FlyBase 5UTR  7529  7679 . +  . gene_symbol "CG11023"; transcript_id "FBtr...
    2L FlyBase exon  7529  8116 . +  . gene_symbol "CG11023"; transcript_id "FBtr...

    Exon Fragment Data for Event Analysis:
    seqname start end   strand gene_id      transcript_id
    2L      7529  8116  +      FBgn0031208  FBtr0300689
    2L      8193  9484  +      FBgn0031208  FBtr0300689
    2L      9839  11344 -      FBgn0002121  FBtr0078169
    """
    all_gtf_columns = ['seqname', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame',
                       'attributes', 'comments']
    drop_columns = ['source', 'feature', 'score', 'frame', 'comments']
    data = pd.read_csv(infile, sep='\t', comment='#', header=None, low_memory=False)
    file_cols = data.columns
    if len(file_cols) < len(all_gtf_columns):
        gtf_cols = all_gtf_columns[:len(file_cols)]
    data.columns = gtf_cols
    drop_cols = [x for x in drop_columns if x in gtf_cols]
    data['seqname'] = data['seqname'].astype(str)
    data['start'] = data['start'].astype(int)
    data['end'] = data['end'].astype(int)
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
    return data


def cli():
    """CLI interface for the 'trand' executable"""
    args = parse_args()
    setup_logging(args.debug, args.verbose, args.log_file)
    logger.debug("Args: {}", args)
    prepare_outdir(args)
    if len(args.infiles) == 1:
        logger.debug("Single file {} analysis", args.ea_mode)
        outfiles = common_outfiles
        if args.ea_mode == "pairwise":
            outfiles.update(pairwise_outfiles)
        else:
            outfiles.update(gene_outfiles)
        try:
            process_single_file(
                args.infiles[0],
                args.ea_mode,
                args.keep_ir,
                args.outdir,
                outfiles,
                args.complexity_only,
                args.skip_plots,
                args.skip_interm,
                args.consolidate,
                args.consol_prefix,
                consol_outfiles
            )
        finally:
            # Only for bedtools. Remove when bedtools are refactored out.
            cleanup()
    else:
        logger.debug("Two files pairwise analysis")
        outfiles = common_outfiles
        outfiles.update(two_gtfs_outfiles)
        outfiles.update(pairwise_outfiles)
        try:
            process_two_files(
                args.infiles,
                args.outdir,
                outfiles,
                args.cpu_cores,
                args.out_pairs,
                args.complexity_only,
                args.skip_plots,
                args.skip_interm,
                args.name1,
                args.name2,
            )
        finally:
            # Only for bedtools. Remove when bedtools are refactored out.
            cleanup()
