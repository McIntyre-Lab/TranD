# -*- coding: utf-8 -*-
"""
Created on Tue Jan 31 12:55:09 2023

@author: k.bankole
"""

"""
Create a customizable upset plot using existing TranD output

Can remove different types of Alternative Splicing
"""

import trand.plot_functions as PF
import trand.io
import argparse
import numpy as np
import pandas as pd


# removing all should just show a grey dot w number of pairs
# default is to show all

def getOptions():
    # Parse command line arguments
    parser = argparse.ArgumentParser(description=
                                     "Make customizable plots using TranD ouput files "
                                     "with the option to ignore specific types of AS."
                                     "Default plot will show " "all types of AS," 
                                     "ignoring all wil show the number of pairs."
                                     "Input a trand output file (csv) (--indir), ignore options "
                                     "(--ignore-3, --ignore-5 --ignoreAD, --ignoreAE, "
                                     "--ignoreIR, --ignoreNSNT), and output path (--outdir)."
                                     )
    

    # Input data
    # Honestly unsure of what the actual input data is... is it a gtf file?
    parser.add_argument(
        "-i",
        "--indir",
        dest="indir",
        required=True,
        help="Location of csv file containing transcriptome AS analysis")
    
    #Ignore Options
    parser.add_argument(
            "-i3",
            "--ignore-3",
            dest="ignore_3",
            action="store_true",
            help="Ignore transcript pairs with 3' variation alternative splicing when creating the Upset Plot"
            )
    
    parser.add_argument(
            "-i5",
            "--ignore-5",
            dest="ignore_5",
            action="store_true",
            help="Ignore transcript pairs with 5' variation alternative splicing when creating the Upset Plot"
            )
    
    parser.add_argument(
            "-iAD",
            "--ignore-AD",
            dest="ignore_AD",
            action="store_true",
            help="Ignore transcript pairs with alternate donor/acceptor alternative splicing when creating the Upset Plot"
            )
    
    parser.add_argument(
            "-iAE",
            "--ignore-AE",
            dest="ignore_AE",
            action="store_true",
            help="Ignore transcript pairs with alternate exon alternative splicing when creating the Upset Plot"
            )
    
    parser.add_argument(
            "-iIR",
            "--ignore-IR",
            dest="ignore_IR",
            action="store_true",
            help="Ignore transcript pairs with intron retention alternative splicing when creating the Upset Plot"
            )
    
    parser.add_argument(
            "-iNSNT",
            "--ignore-NSNT",
            dest="ignore_NSNT",
            action="store_true",
            help="Ignore transcript pairs with no shared nucleotides when creating the Upset Plot"
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
    
    #idea,,, make the prefix the name of the file + whatever was ignored
    # works because the input is already a TranD output file
    # parser.add_argument(
    #     "-x",
    #     "--prefix",
    #     dest="outPrefix",
    #     required=True,
    #     help=(
    #             "Output prefix for plots. "
    #             "Default: no prefix for 1GTF, 'name1_vs_name2' for 2GTF."
    #     )
    # )
    
    parser.add_argument(
        "-f",
        "--force",
        dest="force",
        action="store_true",
        help="Force overwrite existing output directory and files within.",
    )

    args = parser.parse_args()
    return args


# input is a data frame that comes from trand output i think?
# move to plot functions?

#input data frame, legend output directory, ignoreflags
def plot_custom_plot(in_Df, legendOut, ignore_dict):
    pairAS = in_Df[
        [
            "gene_id",
            "transcript_1",
            "transcript_2",
            "flag_alt_exon",
            "flag_alt_donor_acceptor",
            "flag_IR",
            "flag_5_variation",
            "flag_3_variation",
            "flag_no_shared_nt",
            "num_nt_diff",
            "prop_nt_diff",
        ]
    ].copy()
    
    # Ensure number of nt different are int and float values
    pairAS["num_nt_diff"] = pairAS["num_nt_diff"].astype(int)
    pairAS["prop_nt_diff"] = pairAS["prop_nt_diff"].astype(float)
    
    # assures that an exorbitant amt of nucleotides diff are not included in graph
    # bc nsnt will have 100% difference
    pairAS["num_nt_diff"] = np.where(
        pairAS["flag_no_shared_nt"] == 1,
        np.nan,
        pairAS["num_nt_diff"]
    )
    
    # Make AS flags boolean values and rename
    # converts data frame to only AS t/f and number of NT diff and prop diff
    xcrptFlagCols = [
        "transcript_1",
        "transcript_2",
        "flag_alt_exon",
        "flag_alt_donor_acceptor",
        "flag_IR",
        "flag_5_variation",
        "flag_3_variation",
        "flag_no_shared_nt",
    ]
    pairAS[xcrptFlagCols] = pairAS[xcrptFlagCols].astype(bool)
    pairAS = pairAS.rename(
        columns={
            "flag_alt_exon": "Alt. Exon",
            "flag_alt_donor_acceptor": "Alt. Donor/Acceptor",
            "flag_IR": "Intron Retention",
            "flag_5_variation": "5' Variation",
            "flag_3_variation": "3' Variation",
            "flag_no_shared_nt": "No Shared NT",
            "num_nt_diff": "# NT Different",
            "prop_nt_diff": "Proportion\nNT Different",
        }
    )
    
    #set index, sets the index of the dataframe to whatever you want, can have multiple indices
      
    # Create Column Subset based on user arguments
    AScol = []
    for column_header, ignore_flag in ignore_dict.items():
            if not ignore_flag[0]:
                    AScol.append(column_header)
                 
    print(AScol)
    PF.plot_upset(
        pairAS.set_index(AScol),
        "",
        [
            "# NT Different",
            "Proportion\nNT Different",
        ],
    )
    
    # Problems: 
    # 1. What happens when there is only one column, whats the best way to output that?
            # i mean you could just straight print the desired output + num nt diff + prop diff? you dont really
            # need a plot for that
    # 2. What happens when there are no columns? Just output an image of the num nt diff and the prop diff?
    
    # another thing: i dont really understand the legends,,, are they necessary?
        
    

#Builds a dictionary based on user input to customize plot headers    
def ignore_dict_builder(args):
        i_dict = {"3' Variation":[args.ignore_3], "5' Variation":[args.ignore_5], "Alt. Exon":[args.ignore_AE], 
                  "Alt. Donor/Acceptor":[args.ignore_AD], "Intron Retention":[args.ignore_IR], 
                  "No Shared NT":[args.ignore_NSNT]}        
        
        return i_dict
    
    
def main():
        ignore_dictionary = ignore_dict_builder(args)
        
        # trand.io.prepare_outdir(args.outdir, args.force)
        
        inputDf = pd.read_csv(args.indir)
        
        plot_custom_plot(inputDf, args.outdir, ignore_dictionary)
        return 'KSB'


if __name__ == '__main__':
    # Parse command line arguments
    global args
    args = getOptions()
    # ignore_dictionary = ignore_dict_builder(args)
    # testdf = pd.read_csv("/ufgi-vulcana/data/k.bankole/github/TranD/testdata/my_test_csv.csv")
    # print(testdf)
    
    #pairASfix = plot_custom_plot(pd.read_csv(args.indir), args.outdir, ignore_dictionary)

    main()