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
import re
import os
import sys
from importlib.metadata import version
from loguru import logger
from pathlib import Path
from pybedtools import cleanup
from trand.io import prepare_outdir
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
    "ir_fh": "ir_transcripts.csv",
    "ue_fh": "uniq_exons_per_gene.csv",
}
pairwise_outfiles = {"td_fh": "pairwise_transcript_distance.csv"}

gene_outfiles = {
    "er_fh": "event_analysis_er.csv",
    "ef_fh": "event_analysis_ef.csv",
    "ir_fh": "ir_transcripts.csv",
    "ue_fh": "uniq_exons_per_gene.csv",
}
two_gtfs_outfiles = {
    "gtf1_fh": "gtf1_only.gtf",
    "gtf2_fh": "gtf2_only.gtf",
    "md_fh": "minimum_pairwise_transcript_distance.csv",
}


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
        "--subset-pairs",
        action="store",
        dest="subset_file",
        type=str,
        required=False,
        help="""CSV file of transcript pairs that will be processed in pairwise mode and all other pairs will be skipped.
        File should have no header and two columns ordered by transcript_1 and transcript_2. Default is to process all pairs.""",
        default = None
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
        "-x",
        "--prefix",
        action="store",
        type=str,
        required=False,
        help="""Output prefix of various output files. "
        Default: no prefix for 1GTF, 'name1_vs_name2' for 2GTF.""",
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
        "-c",
        "--complexityOnly",
        dest="complexity_only",
        action="store_true",
        help="""Used with 1 or 2 GTF input file(s). Output only transcriptome complexity measures.
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
        input files. Compare pairs of transcripts within a gene. gene - Used with 1 GTF input file.
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
    parser.add_argument(
        "--version",
        action="version",
        version=f"{version('trand')}",
        help="Display version",
    )
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
    regex = re.compile(r"^\w+$", re.ASCII)
    # Validate prefixes
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
        logger.error("Invalid value for the number of CPU cores. Must be 1 or greater.")
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


def setup_logging(debug, verbose, logfile, force):
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
        if force:
            logger.add(logfile, level=level, mode="w")
        else:
            logger.add(logfile, level=level, mode="a")
    logger.debug("Logging level set to : {}", level)


def cli():
    """CLI interface for the 'trand' executable"""
    args = parse_args()
    if not args.outdir:
        args.outdir = str(Path.cwd())
        args.force = True
    setup_logging(args.debug, args.verbose, args.log_file, args.force)
    logger.debug("Args: {}", args)
    prepare_outdir(args.outdir, args.force)
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
                args.cpu_cores,
                args.complexity_only,
                args.skip_plots,
                args.skip_interm,
                args.prefix,
                args.subset_file
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
                args.prefix,
                args.subset_file
            )
        finally:
            # Only for bedtools. Remove when bedtools are refactored out.
            cleanup()
