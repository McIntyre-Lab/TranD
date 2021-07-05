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
import logging
import sys
from loguru import logger
from pathlib import Path
from pybedtools import cleanup
from trand.event_analysis import process_single_file
from trand.event_analysis import process_two_files


# CONFIGURATION
common_outfiles = {'ea_fh': 'event_analysis.csv', 'jc_fh': 'junction_catalog.csv', 'er_fh':
                   'event_analysis_er.csv', 'ef_fh': 'event_analysis_ef.csv'}
pairwise_outfiles = {'td_fh': 'pairwise_transcript_distance.csv'}
gene_outfiles = {'er_fh': 'event_analysis_er.csv', 'ef_fh': 'event_analysis_ef.csv'}
two_gtfs_outfiles = {'gtf1_fh': 'gtf1_only.gtf', 'gtf2_fh': 'gtf2_only.gtf',
                     'md_fh': 'minimum_pairwise_transcript_distance.csv'}


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
        "-c", "--complexityOnly",
        dest='complexity_only',
        action='store_true',
        help="""Output only complexity measures, skipping event analysis and comparison functions
                (default: Perform all analyses and comparisons including complexity calculations)"""
    )
    parser.add_argument(
        "--consolidate",
        dest='consolidate',
        action='store_true',
        help="""Consolidate transcripts with identical junctions prior to evaluation of a single transcriptome
                (remove 5'/3' variation in redundantly spliced transcripts)."""
    )
    parser.add_argument(
        "--consolPrefix",
        dest='consol_prefix',
        default='tr',
        help="""Prefix for consolidated transcript_id values if consolidation is performed (default: tr).
                Prefix must be alphanumeric, can only include \"_\" special character and not contain any spaces."""
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
        "-p", "--pairs",
        type=str,
        choices=['all', 'both', 'first', 'second'],
        dest='out_pairs',
        default='both',
        help="""Output pairwise distance values for:
                all - all transcript pairs in 2 GTF comparison,
                both (default) - only minimum pairs for both datasets,
                first - only minimum pairs for the first dataset,
                second - only minimum pairs for the second dataset"""
    )
    parser.add_argument(
        "-n", "--cpu",
        dest='cpu',
        type=int,
        default=1,
        required=False,
        help="Number of CPUs to use for parallelization (default: 1).",
    )
    parser.add_argument(
        "-1", "--name1",
        dest='name1',
        default="d1",
        required=False,
        help="""For multiple transcriptomes: name of transcriptome for dataset 1, for the first GTF
        file listed, to be used in output (default: \"d1\"). Name must be alphanumeric, can only
        include \"_\" special character and not contain any spaces.""",
    )
    parser.add_argument(
        "-2", "--name2",
        dest='name2',
        default="d2",
        required=False,
        help="""For multiple transcriptomes: name of transcriptome for dataset 2, for the second GTF
        file listed, to be used in output (default: \"d2\"). Name must be alphanumeric, can only
        include \"_\" special character and not contain any spaces.""",
    )
    parser.add_argument(
        "-f", "--force",
        action="store_true",
        help="Force overwrite existing output directory and files within.",
    )
    parser.add_argument(
        "-s", "--skipplots",
        dest='skip_plots',
        action="store_true",
        help="Skip generation of all plots.",
    )
    parser.add_argument(
        "-i", "--skip-intermediate",
        dest='skip_interm',
        action="store_true",
        help="Skip output of intermediate files (such as junction and exon region/fragment files).",
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


def cli():
    """CLI interface for the 'trand' executable"""
    args = parse_args()
    setup_logging(args.debug, args.verbose, args.log_file)
    logger.debug("Args: {}", args)
    infiles = args.infiles
    outdir = handle_outdir(args)
    ea_mode = args.ea_mode
    out_pairs = args.out_pairs
    skip_plots = args.skip_plots
    complexity_only = args.complexity_only
    consolidate = args.consolidate
    if consolidate:
        if [k for k in list(args.consol_prefix) if k.isalnum() or k == "_"] == list(args.consol_prefix):
            consol_prefix = args.consol_prefix
        else:
            logger.error(
                "Invalid consolidated prefix for consolidated transcript_id values: "
                "Must be alphanumerica and can only include '_' special character"
            )
    else:
        consol_prefix = args.consol_prefix
    skip_interm = args.skip_interm
    cpu = args.cpu
    if len(infiles) == 1:
        logger.debug("Single file {} analysis", ea_mode)
        outfiles = common_outfiles
        keep_ir = args.keepir
        if ea_mode == "pairwise":
            outfiles.update(pairwise_outfiles)
        else:
            outfiles.update(gene_outfiles)
        try:
            process_single_file(infiles[0], ea_mode, keep_ir, outdir, outfiles, complexity_only,
                                skip_plots, skip_interm, consolidate, consol_prefix)
        finally:
            cleanup()
    else:
        logger.debug("Two files pairwise analysis")
        outfiles = common_outfiles
        outfiles.update(two_gtfs_outfiles)
        outfiles.update(pairwise_outfiles)
        if [k for k in list(args.name1) if k.isalnum() or k == "_"] == list(args.name1):
            name1 = args.name1
        else:
            logger.error(
                "Invalid name for dataset 1: Must be alphanumeric and can only "
                "include '_' special character"
            )
        if [k for k in list(args.name2) if k.isalnum() or k == "_"] == list(args.name2):
            name2 = args.name2
        else:
            logger.error(
                "Invalid name for dataset 2: Must be alphanumeric and can only "
                "include '_' special character"
            )

        try:
            process_two_files(infiles, outdir, outfiles, cpu, out_pairs, complexity_only,
                              skip_plots, skip_interm, name1, name2)
        finally:
            cleanup()
