#!/usr/bin/env python

import argparse
import pandas as pd
import os
import matplotlib.pyplot as plt
from pathlib import Path

from trand import io
from trand import calculate_complexity as COMP
from trand import plot_functions as PF


def getOptions():
    # Parse command line arguments
    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter,
        description=(
            """
            Calculate complexity measures from gtf (--in-gtf):
                1) Transcripts per gene
                2) Unique exons per gene
                3) Exons per transcript
    
            """
        )
    )

    # Input data
    parser.add_argument(
        "-g",
        "--in-gtf",
        dest="inGTF",
        required=False,
        help=(
                "GTF to calculate complexity plots for."
        )
    )
    
    # Output data
    parser.add_argument(
        "-o",
        "--outdir",
        dest="outdir",
        required=True,
        help=(
                "Output directory, created if missing. "
                "Default: current directory."
        )
    )
    parser.add_argument(
        "-x",
        "--prefix",
        dest="outPrefix",
        required=False,
        help=(
                "Output prefix for plots. "
                "Defaults to GTF file name. "
        )
    )
    parser.add_argument(
        "-f",
        "--force",
        dest="force",
        action="store_true",
        help="Force overwrite existing output directory and files within.",
    )
    
    args = parser.parse_args()
    return args

def main():
    # inGTF = "/nfshome/k.bankole/mnt/exasmb.rc.ufl.edu-blue/mcintyre/share/sex_specific_splicing/fiveSpecies_annotations/fiveSpecies_2_dmel6_ujc.gtf"    
    # outdir = "/nfshome/k.bankole/Desktop"
    # outPrefix = None

    inGTF = args.inGTF
    outdir = args.outdir
    outPrefix = args.outPrefix
    
    if outPrefix is None:
        outPrefix = os.path.basename(inGTF).split('.')[0]
    
    data = io.read_exon_data_from_file(inGTF)
    
    data["num_exon"] = 1
    data["exon_id"] = (
        data["seqname"].map(str)
        + ":"
        + data["start"].map(str)
        + ":"
        + data["end"].map(str)
        + ":"
        + data["strand"].map(str)
    )
    
    transcriptDF = (
        data.groupby("transcript_id").agg({"num_exon": "sum"}).reset_index()
    )
    
    geneDF = (
        data.groupby("gene_id")
        .agg({"transcript_id": "nunique", "num_exon": "sum", "exon_id": "nunique"})
        .reset_index()
        .rename(columns={"transcript_id": "num_transcript", "exon_id": "num_uniq_exon"})
    )
    
    PF.plot_complexity_box(
        geneDF,
        transcriptDF,
        outdir,
        "{}/{}_complexity_plots.rtf".format(outdir, outPrefix),
    )
    plt.savefig(
        "{}/{}_complexity_plots.png".format(outdir, outPrefix),
        dpi=600,
        format="png",
    )
    plt.clf()
        
    
    

if __name__ == '__main__':
    # Parse command line arguments
    global args
    args = getOptions()
    main()
