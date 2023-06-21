#!/usr/bin/env python

import argparse
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

# Import trand plotting functions
from trand import plot_two_gtf_pairwise as P2GP
from trand import plot_one_gtf_pairwise as P1GP
from trand import plot_one_gtf_gene as P1GG

# Import trand io functions
from trand import io

def getOptions():
    # Parse command line arguments
    parser = argparse.ArgumentParser(
        description=(
                "Make TranD plots from output files for 1 GTF gene "
                "(--exon-region, --exon-EF, --intron-retention, "
                "--uniq-exon), 1 GTF pairwise (--pairwise-1gtf), or "
                "2 GTF pairwise (--pairwise-2gtf, --name1, --name2, --gtf1, "
                "--gtf2). If both 1 GTF gene exon EF file "
                "(--exon-EF) and 1 GTF pairwise distance file "
                "(--pairwise-1gtf) are included, a new density plot of "
                "nucleotide variablility separated by AS type will be made."
        )
    )

    # Input data
    parser.add_argument(
        "-er",
        "--exon-region",
        dest="er_data",
        required=False,
        help=(
                "For 1 GTF gene plots: "
                "Exon regions output file (event_analysis_er.csv) from "
                "TranD (must be used with --exon-EF argument)."
        )
    )
    parser.add_argument(
        "-ef",
        "--exon-EF",
        dest="ef_data",
        required=False,
        help=(
                "For 1 GTF gene plots: "
                "Exon EF output file (event_analysis_ef.csv) from "
                "TranD (must be used with --exon-region argument)."
        )
    )
    parser.add_argument(
        "-ir",
        "--intron-retention",
        dest="ir_data",
        required=False,
        help=(
                "For 1 GTF gene plots: "
                "Intron retention (IR) transcripts output file "
                "(ir_transcripts.csv)."
        )
    )
    parser.add_argument(
        "-ue",
        "--uniq-exon",
        dest="ue_data",
        required=False,
        help=(
                "For 1 GTF gene plots: "
                "Unique exons per gene output file "
                "(uniq_exons_per_gene.csv)."
        )
    )
    parser.add_argument(
        "-p1",
        "--pairwise-1gtf",
        dest="td_data1",
        required=False,
        help=(
                "For 1 GTF pairwise plots: "
                "Transcript distance output file "
                "(pairwise_transcript_distance.csv) from TranD 1 GTF pairwise."
        )
    )
    parser.add_argument(
        "-t",
        "--density-threshold",
        dest="density_th",
        required=False,
        type=float,
        help=(
                "For 1 GTF additional plot (combining AS and nucleotide "
                "variability): KDE density threshold value to use as the "
                "maximum for the Y-axis of the gene density plot."
        )
    )
    parser.add_argument(
        "-p2",
        "--pairwise-2gtf",
        dest="td_data2",
        required=False,
        help=(
                "For 2 GTF pairwise plots: "
                "Transcript distance output file "
                "(pairwise_transcript_distance.csv or "
                "minimum_pairwise_transcript_distance.csv) from TranD 2 "
                "GTF pairwise."
        )
    )
    parser.add_argument(
        "-n1",
        "--name1",
        dest="name1",
        required=False,
        default="d1",
        help=(
                "For 2 GTF pairwise distance file, name used for dataset 1 "
                "that is appended to transcript_1 values "
                "(default = d1)."
        )
    )
    parser.add_argument(
        "-n2",
        "--name2",
        dest="name2",
        required=False,
        default="d2",
        help=(
                "For 2 GTF pairwise distance file, name used for dataset 2 "
                "that is appended to transcript_2 values "
                "(default = d2)."
        )
    )
    parser.add_argument(
        "-g1",
        "--gtf1",
        dest="gtf1",
        required=False,
        help=(
                "For 2 GTF pairwise plots: "
                "TranD output GTF file for genes only "
                "in dataset 1 (gtf1_only.gtf). While not required, "
                "with no GTF provided, pie charts of gene counts will "
                "show 0 genes exclusive to one dataset."
        )
    )
    parser.add_argument(
        "-g2",
        "--gtf2",
        dest="gtf2",
        required=False,
        help=(
                "For 2 GTF pairwise plots: "
                "TranD output GTF file for genes only "
                "in dataset 2 (gtf2_only.gtf). While not required, "
                "with no GTF provided, pie charts of gene counts will "
                "show 0 genes exclusive to one dataset."
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
        required=True,
        help=(
                "Output prefix for plots. "
                "Default: no prefix for 1GTF, 'name1_vs_name2' for 2GTF."
        )
    )
    parser.add_argument(
        "-f",
        "--force",
        dest="force",
        action="store_true",
        help="Force overwrite existing output directory and files within.",
    )
    
    parser.add_argument(
            "-ns",
            "--no-suffix",
            dest="noSuff",
            action="store_true",
            help="Add this argument if there is no suffix on the transcripts in the 2GTF data"
            )
    
    args = parser.parse_args()
    return args


er_df_cols = ['gene_id', 'er_id', 'er_chr', 'er_start', 'er_end', 'er_strand', 'er_exon_ids',
              'er_transcript_ids', 'gene_transcript_ids', 'exons_per_er', 'transcripts_per_er',
              'transcripts_per_gene',  'er_ir_flag', 'er_annotation_frequency']
ef_df_cols = ['gene_id', 'er_id', 'ef_id', 'ef_chr', 'ef_start', 'ef_end', 'ef_strand',
              'ef_exon_ids', 'ef_transcript_ids', 'exons_per_ef', 'transcripts_per_ef',
              'ef_ir_flag', 'ea_annotation_frequency']
td_df_cols = ['gene_id','transcript_1','transcript_2','num_jxn_only_T1','num_jxn_only_T2',
             'num_jxn_ovlp','prop_jxn_noOvlp','prop_jxn_ovlp',
             'jxn_only_T1','jxn_only_T2','jxn_same','num_ER_only_T1',
             'num_ER_only_T2','num_ER_ovlp','prop_ER_noOvlp','prop_ER_ovlp',
             'ER_only_T1','ER_only_T2','ER_ovlp','num_EF_only_T1',
             'num_EF_only_T2','num_EF_ovlp','prop_EF_noOvlp',
             'prop_EF_ovlp','EF_only_T1','EF_only_T2','EF_ovlp',
             'num_exon_only_T1','num_exon_only_T2',
             'num_exon_ovlp','num_IR_EF_T1','num_IR_EF_T2',
             'IR_EF_T1','IR_EF_T2','num_nt_ovlp','num_nt_only_T1',
             'num_nt_only_T2','num_nt_noOvlp','total_nt','prop_nt_noOvlp','prop_nt_ovlp','num_nt_only_T1_in_ovlpER',
             'num_nt_only_T2_in_ovlpER','num_nt_ovlp_in_ovlpER','total_nt_in_ovlpER',
             'prop_nt_noOvlp_in_ovlpER','prop_nt_ovlp_in_ovlpER','num_nt_only_T1_in_uniqER',
             'num_nt_only_T2_in_uniqER','flag_FSM','flag_IR','flag_5_var',
             'flag_3_var','flag_alt_DA','flag_alt_exon','flag_no_ovlp_nt']
ir_df_cols = ['er_transcript_ids']
ue_df_cols = ['gene_id', 'num_uniq_exon']

def check_args(args):
    """
    Check input arguments to determine which plots can be made
    """
    ntVar = False
    oneGTFpair = False
    oneGTFgene = False
    twoGTF = False
    # Check for 1 GTF pairwise outupt and 1 GTF gene exon EF output
    #   for nt variablility density plot
    if args.td_data1 is not None and args.ef_data is not None:
        ntVar = True
    # Check for 1 GTF pairwise TranD output
    if args.td_data1 is not None:
        oneGTFpair = True
    # Check for 1 GTF gene TranD output
    if (
        (args.er_data is not None) and 
        (args.ef_data is not None) and
        (args.ir_data is not None) and
        (args.ue_data is not None)        
        ):
        oneGTFgene = True
    elif (
        (args.er_data is not None) or 
        (args.ef_data is not None) or
        (args.ir_data is not None) or
        (args.ue_data is not None)
        ) and (args.td_data1 is None):
        exit(
                "!!!ERROR: Must provide exon region, exon EF, "
                "intron retention, and unique exons per gene files to "
                "properly plot TranD 1 GTF gene plots."
        )
    if args.td_data2 is not None:
        twoGTF = True
        if args.name1 is None:
            print("!!!WARNING: Using default name for dataset 1 - \"d1\".")
        if args.name2 is None:
            print("!!!WARNING: Using default name for dataset 2 - \"d2\".")
        if args.gtf1 is None:
            print(
                    "!!!WARNING: No GTF file provided for genes exclusive to "
                    "dataset 1 (gtf1_only.gtf). Plots will not reflect "
                    "accurate counts of genes exclusive to one dataset."
            )
        if args.gtf2 is None:
            print(
                    "!!!WARNING: No GTF file provided for genes exclusive to "
                    "dataset 2 (gtf2_only.gtf). Plots will not reflect "
                    "accurate counts of genes exclusive to one dataset."
            )
    return ntVar, oneGTFpair, oneGTFgene, twoGTF


def validate_input(inType, args):
    if inType == "oneGTFpair" or inType == "ntVar":
        td_data = pd.read_csv(args.td_data1, low_memory=False)
        # Check that all transcript distance columns are present
        if [c for c in td_data if c in td_df_cols] != td_df_cols:
            exit(
                    "Single GTF pairwise transcript distance file missing "
                    "columns from TranD output."
            )
        else:
            if inType == "ntVar":
                # Check if all exon EF columns are present
                ef_data = pd.read_csv(args.ef_data, low_memory=False)
                if [c for c in ef_data if c in ef_df_cols] == ef_df_cols:
                    return td_data, ef_data
                else:
                    exit(
                        "!!!ERROR: Single GTF exon EF file missing "
                        "columns from TranD output."
                    )
            else:
                return td_data
    elif inType == "oneGTFgene":
        er_data = pd.read_csv(args.er_data, low_memory=False)
        ef_data = pd.read_csv(args.ef_data, low_memory=False)
        ir_data = pd.read_csv(args.ir_data, low_memory=False)
        ue_data = pd.read_csv(args.ue_data, low_memory=False)
        # Check that all exon region and exon EF columns are present
        if (
                [c for c in er_data if c in er_df_cols] == er_df_cols and
                [c for c in ef_data if c in ef_df_cols] == ef_df_cols and
                [c for c in ir_data if c in ir_df_cols] == ir_df_cols and
                [c for c in ue_data if c in ue_df_cols] == ue_df_cols
            ):
            return er_data, ef_data, ir_data, ue_data
        else:
            if [c for c in er_data if c in er_df_cols] != er_df_cols:
                exit(
                        "!!!ERROR: Single GTF exon region file missing "
                        "columns from TranD output."
                )
            if [c for c in ef_data if c in ef_df_cols] != ef_df_cols:
                exit(
                        "!!!ERROR: Single GTF exon EF file missing "
                        "columns from TranD output."
                )
            if [c for c in ir_data if c in ir_df_cols] != ir_df_cols:
                exit(
                        "!!!ERROR: Single GTF intron retention file missing "
                        "columns from TranD output."
                )
            if [c for c in ue_data if c in ue_df_cols] != ue_df_cols:
                exit(
                        "!!!ERROR: Single GTF unique exons per gene file "
                        "missing columns from TranD output."
                )
    elif inType == "twoGTF":
        td_data = pd.read_csv(args.td_data2, low_memory=False)
        if args.gtf1 is not None:
            try:
                f1_odds = io.read_exon_data_from_file(args.gtf1)
            except:
                f1_odds = pd.DataFrame(columns=['gene_id','transcript_id'])
                print("!!!WARNING: GTF with genes exclusive to dataset 1 is empty.")
        else:
            f1_odds = pd.DataFrame(columns=['gene_id','transcript_id'])
        if args.gtf2 is not None:
            try:
                f2_odds = io.read_exon_data_from_file(args.gtf2)
            except:
                f2_odds = pd.DataFrame(columns=['gene_id','transcript_id']) 
                print("!!!WARNING: GTF with genes exclusive to dataset 2 is empty.")
        else:
            f2_odds = pd.DataFrame(columns=['gene_id','transcript_id'])
        name1 = args.name1
        name2 = args.name2
        # Check that all transcript distance columns are present
        if [c for c in td_data if c in td_df_cols] == td_df_cols:
            # Check that name1 and name2 provided are present
            #   in the transcript_1 and transcript 2 columns
            if not args.noSuff: 
                    if td_data['transcript_1'].str.contains(args.name1).all():
                            if td_data['transcript_2'].str.contains(args.name2).all():
                                    return td_data, name1, name2, f1_odds, f2_odds
            else:
                    return td_data, name1, name2, f1_odds, f2_odds

        else:
            exit(
                    "!!!ERROR: Two GTF pairwise transcript distance file "
                    "missing columns from TranD output."
            )
    else:
        exit("Invalid input type.")
def split_column_by_sep(df, col_name=None, sep=None, sort_list=None):
    """
    Split variable by some character like '|' or ',' and
    keep all other values the same
    """
    if col_name == None:
        col_name = 'transcript_id'
    if sep == None:
        sep = "|"
    splitList = df[col_name].str.split(sep).apply(pd.Series, 1).stack()
    splitList.index = splitList.index.droplevel(-1)
    tempDF = df.copy()
    del(tempDF[col_name])
    splitDF = tempDF.join(splitList.rename(col_name))
    if sort_list != None:
        splitDF = splitDF.sort_values(by=sort_list)
    del(tempDF, splitList)
    return splitDF

def plot_gene_prop_nt_variablility_AS(td_data, ef_data, density_th=None):
    """
    Plot distribution of proportion of nucleotide variability across all
    multi-transcript genes, split by AS categories
    """
    # Get lengths of exon EFs
    ef_data['num_nt_varying'] = np.where(
            ef_data['ea_annotation_frequency']!="constitutive",
            ef_data['ef_end'] - ef_data['ef_start'],
            0,
    )
    ef_data['num_nt_total'] = ef_data['ef_end'] - ef_data['ef_start']
    # For each gene, get the number of varying nucleotides (not constitutive)
    ef_gene_data = ef_data.groupby('gene_id')[
            [
                    'num_nt_varying',
                    'num_nt_total'
            ]
    ].sum().reset_index()
    # Get propotion of varying nucleotides
    ef_gene_data['prop_varying_nt'] = (
            ef_gene_data['num_nt_varying'] / ef_gene_data['num_nt_total']
    )
    # Get only multi-transcript genes
    ef_data_split = split_column_by_sep(
            ef_data,
            col_name="ef_transcript_ids",
            sep="|",
    )
    xcrpt_per_gene = ef_data_split.groupby('gene_id')[
            'ef_transcript_ids'
    ].nunique().reset_index()
    multi_xcrpt_gene = xcrpt_per_gene[
            xcrpt_per_gene['ef_transcript_ids']>1
    ]
    ef_gene_data = ef_gene_data[
            ef_gene_data['gene_id'].isin(
                    multi_xcrpt_gene['gene_id']
            )
    ]

    # Plot densities of genes with each category of AS (not mutually exclusive)
    td_gene_data = td_data.groupby('gene_id')[
            [c for c in td_data.columns if "flag_" in c]
    ].max().reset_index()
    # Check if IR were dropped from EF file (no -k option used in trand)
    if ef_data['ef_ir_flag'].sum() == 0:
        labelDict = {
                'flag_5_var': "5\' Variation",
                'flag_3_var': "3\' Variation",
                'flag_alt_DA': "Alt. Donor/Acceptor",
                'flag_alt_exon': "Alt. Exon"
        }
    else:
        labelDict = {
                'flag_IR': "Intron Retention",
                'flag_5_var': "5\' Variation",
                'flag_3_var': "3\' Variation",
                'flag_alt_DA': "Alt. Donor/Acceptor",
                'flag_alt_exon': "Alt. Exon"
        }
    colorDict = {
            'flag_IR': '#1f77b4',
            'flag_5_var': '#ff750e',
            'flag_3_var': '#2ca02c',
            'flag_alt_DA': '#d62728',
            'flag_alt_exon': '#9467bd'
    }
    for col in labelDict.keys():
        genes = td_gene_data[td_gene_data[col]==1]['gene_id']
        sns.distplot(
#                ef_gene_data[ef_gene_data['gene_id'].isin(genes)]['prop_varying_nt'],
                ef_gene_data[
                        (ef_gene_data['gene_id'].isin(genes))&
                        (ef_gene_data['prop_varying_nt']!=0)
                ]['prop_varying_nt'],
                hist=False,
                label=labelDict[col]+", n = ({})".format(
                        len(ef_gene_data[ef_gene_data['gene_id'].isin(genes)])
                        ),
                color=colorDict[col],
        )
    plt.xlim(0,1)
    if density_th is not None:
        plt.ylim(0, density_th)
    plt.ylabel("Density")
    plt.xlabel("Proportion of Variable Nucleotides")
    plt.legend()
    plt.tight_layout()

def main():
    # Check input arguments
    ntVar, oneGTFpair, oneGTFgene, twoGTF = check_args(args)

    # Prepare outdir
    io.prepare_outdir(args.outdir, args.force)
    outdir = args.outdir

    # Validate inputs and plot
    if ntVar:
        td_data, ef_data = validate_input(inType="ntVar", args=args)
        plot_gene_prop_nt_variablility_AS(td_data, ef_data, args.density_th)
        plt.savefig(
                "{}/gene_AS_prop_nt_variablility.png".format(outdir),
                dpi=600,
                format="png"
        )
        plt.clf()
    if oneGTFpair:
        td_data = validate_input(inType="oneGTFpair", args=args)
        P1GP.plot_one_gtf_pairwise(
                outdir,
                td_data,
                prefix=args.outPrefix,
        )
    if oneGTFgene:
        er_data, ef_data, ir_data, ue_data = (
            validate_input(inType="oneGTFgene", args=args)
        )
        P1GG.plot_one_gtf_gene(
                er_data,
                ef_data,
                ir_data,
                ue_data,
                outdir,
                prefix=args.outPrefix,
        )
    if twoGTF:
            try:
                    td_data, name1, name2, f1_odds, f2_odds = (
                            validate_input(inType="twoGTF", args=args)
                    )
            except TypeError:
                    raise TypeError("Input was not validated. You may be missing a -ns in your arguments.")
            if args.outPrefix is None:
                    args.outPrefix = str(name1) + "_vs_" + str(name2)
                
            P2GP.plot_two_gtf_pairwise(
                outdir,
                td_data,
                f1_odds,
                f2_odds,
                name1=name1,
                name2=name2,
                prefix=args.outPrefix
           )


if __name__ == '__main__':
    # Parse command line arguments
    global args
    args = getOptions()
    main()
