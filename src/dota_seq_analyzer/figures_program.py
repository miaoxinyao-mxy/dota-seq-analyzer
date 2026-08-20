# imports
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.colors as colors
import pandas as pd
from typing import List, Dict, Tuple
from collections import Counter
from helper_functions import get_arg_names, ensure_output_directories
import os
import argparse

# ===============================================================================================

# Primary function for making all 3 figures

# 2026-08-10: Use 30 cells as the default minimum for ASV-level visualization.
# Reason: a 40-cell cutoff removes all ASV groups from moderate-sized single-cell datasets.
def make_figures(
    use_asvs_str: str, unfiltered_barcode_summary_tsv: str, final_asv_barcode_summary_tsv: str, 
    asv_barcode_summary_no_sub_args_tsv: str, primers_file: str, 
    b_with_ids: str, asv_arg_table_tsv: str,  
    asv_arg_figure: str, barcode_group_size_figure: str, primer_balance_figure: str,
    first_gene_column_num: int, global_mle_tax_tsv: str, global_asv_tsv: str = None,
    arg_threshold: float = 0.01, min_cells_per_asv: int = 30, figure_dpi: int = 300):

    use_asvs = determine_use_asvs(use_asvs_str) # convert use_asvs from yes/no into a bool True/False value

    # 0) Non-figure related: write global taxonomic classifications file
    write_global_tax_classification_file(final_asv_barcode_summary_tsv, global_mle_tax_tsv)

    # 1) ASV-ARG Table
    make_asv_arg_table(asv_arg_figure, final_asv_barcode_summary_tsv, asv_arg_table_tsv, 
        min_cells_per_asv, arg_threshold, figure_dpi, first_gene_column_num, 
        use_asvs, global_mle_tax_tsv, global_asv_tsv)
    
    # 2) Barcode Group Size QC
    make_barcode_group_size_figure(barcode_group_size_figure, 
        unfiltered_barcode_summary_tsv, primers_file, b_with_ids, figure_dpi)
    
    # 3) Primer Balance (ARG vs 16s primers) QC 
    make_primer_balance_figure(primer_balance_figure, asv_barcode_summary_no_sub_args_tsv, primers_file, figure_dpi)

# ===============================================================================================

# Core functions for making figures

def make_asv_arg_table(
    asv_arg_figure: str, final_asv_barcode_summary_tsv: str, asv_arg_table_tsv: str, 
    min_cells_per_asv: int, arg_threshold: float, figure_dpi: int, first_gene_column_num: int, 
    use_asvs: bool, global_mle_tax_tsv, global_asv_tsv: str = None):
    """
    Figure: Heat map of ASVs with ARGs
    Purpose: determine which ARGs each "species" has (where each "species" is identified with its ASV)
    """
    # 2026-08-10: Present the existing ASV association plot as a generic target plot.
    # Reason: target panels may contain AMR, marker genes, or phase-variation loci.
    df_asv_arg, final_gene_names = create_asv_arg_matrix(
        use_asvs, final_asv_barcode_summary_tsv,  
        asv_arg_table_tsv, min_cells_per_asv, first_gene_column_num, 
        global_mle_tax_tsv, global_asv_tsv)

    if use_asvs:
        y_axis_label = "ASV"
        # set labels equal to list of asv names, with taxonomic classifications added to each asv label
        asv_or_tax_labels = []
        for idx in df_asv_arg.index.to_list():
            most_specific_tax = ""
            full_tax = df_asv_arg.loc[idx, "Predicted taxonomy"].split(" | ")
            for tax_lvl in full_tax:
                if tax_lvl.split(" - ")[1] != "None":
                    most_specific_tax = tax_lvl
            asv_or_tax_labels.append(f"{idx}: {most_specific_tax}") 

        # remove taxonomy column, now that we're done using it
        df_asv_arg.drop("Predicted taxonomy", axis = 1, inplace = True) 

    else:
        y_axis_label = "MLE Taxonomic Classification"
        # set labels equal to list of tax names, edited to include only the most specific taxonomy level
        asv_or_tax_labels = []
        for idx in df_asv_arg.index.to_list():
            most_specific_tax = ""
            full_tax = idx.split(" | ")
            for tax_lvl in full_tax:
                if tax_lvl.split(" - ")[1] != "None":
                    most_specific_tax = tax_lvl
            asv_or_tax_labels.append(most_specific_tax)

    # adjust length of long arg name to better fit the x axis
    for i in range(len(final_gene_names)):
        if final_gene_names[i] == "OXA-48 like (oxa-3)":
            final_gene_names[i] = "OXA-48 like\n(oxa-3)"

    # edit the ASV-ARG table - nullify values below arg threshold
    df_asv_arg.mask(df_asv_arg < arg_threshold, 0, inplace = True)

    # 2026-08-10: Produce an explanatory figure when no heatmap values remain.
    # Reason: small datasets or strict thresholds can yield an empty/all-zero matrix, which Matplotlib cannot plot with LogNorm.
    if df_asv_arg.empty or not np.any(df_asv_arg.to_numpy() > 0):
        fig, axs = plt.subplots(figsize=(8, 3))
        axs.axis("off")
        axs.text(
            0.5,
            0.5,
            "No taxonomic groups met the cell and target thresholds.",
            ha="center",
            va="center",
        )
        fig.savefig(asv_arg_figure, bbox_inches="tight", pad_inches=0.3, dpi=figure_dpi)
        plt.close(fig)
        return

    # plot the ASV-ARG table
    axs = plt.matshow(df_asv_arg.to_numpy(), norm=colors.LogNorm(), cmap = "Blues").axes
    axs.set_xlabel("Target", fontweight = "bold", fontsize = 11)
    axs.xaxis.tick_bottom()
    axs.set_ylabel(y_axis_label, fontweight = "bold", fontsize = 11)
    axs.set_xticks(np.arange(0, len(final_gene_names), 1), final_gene_names, rotation = "vertical") 
    axs.set_yticks(np.arange(0, len(asv_or_tax_labels), 1), asv_or_tax_labels)
    axs.set_title("Targets Present in the Most Common Taxonomic Classifications", fontweight = "bold", fontsize = 15, pad = 30)
    plt.colorbar()
    plt.savefig(asv_arg_figure, bbox_inches = "tight", pad_inches = 0.3, dpi = figure_dpi)
    plt.show()

def make_barcode_group_size_figure(
    barcode_group_size_figure, unfiltered_barcode_summary_tsv, 
    primers_file, b_with_ids, figure_dpi):
    # Figure: barcode group size; i.e. # of reads for barcodes
    # Purpose: filtering; determine threshold min # of reads to use when filtering out barcodes, 
    # such that most of the dataset's reads are still preserved
    # Adapted code from original dota-seq paper


    # obtain list of sizes of all barcode groups (i.e. # of reads per barcode)
    # there will be three graphs: 1) only 16s reads, 2) all classified reads (16s + ARG), 3) all reads (16s + ARG + unclassified)
    sum_16s_reads = get_num_16s_reads_for_barcodes(unfiltered_barcode_summary_tsv)
    sum_classified_reads = get_num_reads_for_classified_barcodes(unfiltered_barcode_summary_tsv, primers_file)
    sum_all_reads = get_num_reads_for_all_barcodes(b_with_ids)

    # plot the graphs

    fig, axs = plt.subplots(1, 3, figsize = [9, 5])
    fig.tight_layout()    

    ax = axs[0]
    jackpottocurve(fig, ax, sum_16s_reads, "Only 16s", vline=5)
    ax = axs[1]
    jackpottocurve(fig, ax, sum_classified_reads, "Classified (16S + targets)", vline=0)
    ax = axs[2]
    jackpottocurve(fig, ax, sum_all_reads, "All (16S + targets + unclassified)", vline=0)

    plt.savefig(barcode_group_size_figure, bbox_inches = "tight", pad_inches = 0.3, dpi = figure_dpi)

def make_primer_balance_figure(
    primer_balance_figure, asv_barcode_summary_no_sub_args_tsv, primers_file, figure_dpi):
    # Figure: ARG vs 16s primer balance
    # Purpose: determine if the ratio of ARG to 16s primers needs to be adjusted, for any of the ARGs
    # Adapted code from original dota-seq paper

    c_palette = ["#EE7032", "#FC8609", "#F5A232", "#F7DD48", "#C4E54C","#8CC860",
            "#6CC860", "#5ACCAA", "#52CECE", "#4DAEEE", "#5877D6", "#C186F1", "#D441D6"]

    # ratio calculation: (# of ARG#1 reads) / (# of ARG#1 reads + # of 16s reads)  
    # (where "ARG#1" represents whichever of the 23 ARGs you're currently looking at)
    # note that this ratio calculation is done only for the barcodes that are classified as having that specific ARG
    arg_ratios = get_ARG_to_16s_ratios(asv_barcode_summary_no_sub_args_tsv, primers_file)

    # plot the graphs
    fig, axs = plt.subplots(len(arg_ratios), 1, figsize = (10, len(arg_ratios)))
    fig.tight_layout()
    i = 0
    for arg in arg_ratios:
        ax = axs[i]
        ax.hist(arg_ratios[arg], bins = 100, range = (0.1, 1), color = c_palette[i % 13])
        ax.set_yscale("log")
        ax.set_title(arg)
        ax.set_xticks([0,1], ["16S-only", "target-only"])
        #ax.set_yticklabels([0])
        i += 1

    plt.savefig(primer_balance_figure, bbox_inches = "tight", pad_inches = 0.3, dpi = figure_dpi)

# ===============================================================================================

# Non-figure related

# Write global taxonomic classifications file
# Done now as opposed to earlier when making barcode summary, to account for barcodes that have been filtered out

def write_global_tax_classification_file(
    barcode_summary_tsv: str, global_mle_tax_tsv: str): 
    """
    Extract all unique MLE taxonomic classifications from the current barcode summary file.
    Write each taxonomic classification to a file, along with its corresponding cell count,
     in order of most to least common taxonomic classifications.
    Create this file now, as opposed to earlier when making the barcode summary. Reason: some barcodes 
     have been filtered out during the intermediate steps, which will have changed the cell count 
     for many taxonomic classifications. So, now is an ideal time to make this file.
    Purpose: this file will be used to get the most common taxonomic classifications, for the heat map figure
     (used when user chooses not to do ASV typing).
    """
    df = pd.read_csv(barcode_summary_tsv, sep="\t", index_col = "Barcode")
    global_tax_counter = Counter(df["Predicted taxonomy"].to_list())
    with open(global_mle_tax_tsv, 'w') as f:
        f.write("MLE_Taxonomic_Classification\tCell_Count\n")
        for tax_classification, cell_count in global_tax_counter.most_common():
            f.write(f"{tax_classification}\t{cell_count}\n")


# Determine whether or not to use ASVs -> convert yes/no into True/False
def determine_use_asvs(use_asvs_str: str):
    assert use_asvs_str in ["yes", "no"]
    if use_asvs_str == "yes":
        return True
    elif use_asvs_str == "no":
        return False

# ===============================================================================================

# Helper functions for ASV-ARG Table Figure

def get_asv_names(
    global_asv_file: str, min_cells_per_asv: int = 10) -> List[str]:
    """Extract the names of all ASVs that have the desired minimum # of cells associated with them"""
    asv_names = []
    with open(global_asv_file, 'r') as f:
        line = f.readline() # header row
        line = f.readline() # first content row
        while line != "":
            name, num_cells, seq = line.strip().split("\t")
            if int(num_cells) >= min_cells_per_asv:
                asv_names.append(name)
                line = f.readline()
            else:
                break # assume global_asv is sorted from most to least common asvs
                      # so if current asv has < min_cells_per_asv,
                      # the rest of the asvs will also be below threshold, so break the loop
    return asv_names

def get_tax_names(
    global_mle_tax_tsv: str, min_cells_per_tax_classification: int) -> List[str]:
    """Extract all taxonomic classifications that have the desired minimum # of cells associated with them"""
    tax_names = []
    with open(global_mle_tax_tsv, 'r') as f:
        line = f.readline() # header row
        line = f.readline() # first content row
        while line != "":
            name, num_cells = line.strip().split("\t")
            if int(num_cells) >= min_cells_per_tax_classification:
                tax_names.append(name)
                line = f.readline()
            else:
                break # assume global_mle_tax_file is sorted from most to least common taxonomic classifications
                      # so if current tax_classification has < min_cells_per_tax_classification,
                      # the rest of the tax classifications will also be below threshold, so break the loop
    return tax_names

def summarize_ASV_with_ARG(
    filtered_counts_summary_tsv: str, asv_or_tax_names: List[str], 
    first_gene_column_num: int, use_asvs: bool) -> Tuple:
    """
    Create and return a dataframe relating either ASVs or MLE taxonomic classifications, with their present ARGs.
    Specifically, record the fraction of cells associated with each ASV/tax, that have each of the given ARGs.
    Note that the input barcode summary is already filtered by ARGs, and if applicable, filtered by sub-ARGs as well.
    This removes background noise when determining which ASVs/taxes have which ARG.
    """

    # load ASV barcode summary, and extract all the final ARG & sub-ARG names from the barcode summary
    # note that gene columns start at "first_gene_column_num" (where this # is based on zero-indexed columns, 
    # and given that the "Barcode" column is not considered a column, but instead just the df index)
    # and these gene columns go until the last column of the barcode summary

    df_summary = pd.read_csv(filtered_counts_summary_tsv, sep="\t", index_col = "Barcode")
    all_col_names = df_summary.columns.to_list()
    final_gene_names = all_col_names[first_gene_column_num:] 
    
    if use_asvs:
        df_asv_arg = pd.DataFrame(columns = (["Predicted taxonomy"] + final_gene_names))
        col = "Assigned_core_asv"
    else:
        df_asv_arg = pd.DataFrame(columns = final_gene_names)
        col = "Predicted taxonomy"

    for asv_or_tax in asv_or_tax_names:
        df_one_asv_or_tax = df_summary[df_summary[col] == asv_or_tax]

        # determine ARG fractions
        num_cells = len(df_one_asv_or_tax)
        arg_fractions_row = [] # record the fraction of cells for the given ASV that have each of the given ARGs
                               # e.g. [0.08, 0.37, 0.05, ...] if ASV #1 had 8% of its cells having ARG#1, 
                               # 37% of its cells with ARG#2, 5% of its cells with ARG#3, etc.
        for arg in final_gene_names:
            arg_fractions_row.append(len(df_one_asv_or_tax[df_one_asv_or_tax[arg] != 0]) / num_cells)
        
        if use_asvs:
            # determine MLE taxonomy corresponding to that ASV
            tax_counts = Counter(df_one_asv_or_tax["Predicted taxonomy"].to_list())
            taxonomy = tax_counts.most_common(1)[0][0]      
            row = [taxonomy] + arg_fractions_row
        else:
            row = list(arg_fractions_row)

        df_asv_arg.loc[asv_or_tax] = row

    return df_asv_arg, final_gene_names
        
def create_asv_arg_matrix(
    use_asvs: bool, filtered_counts_summary_tsv: str, 
    asv_arg_tsv: str, min_cells_per_asv_or_tax: int, first_gene_column_num: int,
    global_mle_tax_tsv: str, global_asv_file: str = None):
    """
    Main function to create the ASV-ARG table
    Creates and returns a dataframe relating ASVs with their present ARGs
    Also writes this dataframe to a tsv file.
    Note that this is called "asv-arg matrix", but it can be used without ASVs as well;
    in that case, MLE taxonomic classifications would be used as a replacement for ASVs.
    """
    if use_asvs:
        asv_names = get_asv_names(global_asv_file, min_cells_per_asv_or_tax)
        df_asv_arg, final_gene_names = summarize_ASV_with_ARG(filtered_counts_summary_tsv, 
                                       asv_names, first_gene_column_num, use_asvs)
        label = "ASV"
    else:
        tax_names = get_tax_names(global_mle_tax_tsv, min_cells_per_asv_or_tax)
        df_asv_arg, final_gene_names = summarize_ASV_with_ARG(filtered_counts_summary_tsv, 
                                       tax_names, first_gene_column_num, use_asvs)
        label = "MLE Taxonomic Classification"

    df_asv_arg.to_csv(asv_arg_tsv, sep = "\t", index_label = label)

    return df_asv_arg, final_gene_names

# ===============================================================================================

# Helper functions for Barcode Group Size Figure

def jackpottocurve(fig, ax, barsize, type_of_reads, vline):
    #jackpottocurve is a cumulative histogram of the read groups based on how many reads are in the group

    #prepare data to plot the jackpottocurve
    X = np.array(barsize) #make an np array of the data so we can use the various np functions
    X.sort()
    X_lorenz = X.cumsum() / X.sum() #turns it into a cumsum percentage (for y-axis)
    X_lorenz = np.insert(X_lorenz, 0, 0) #insert the 0,0 point
    X = np.insert(X, 0, 0) #insert 0,0 point for X as well

    # scatter plot of the jackpottocurve (modification of the lorenz curve)
    ax.scatter(X, X_lorenz, 
                marker='o', color='indigo', s=100) #x axis is cumulative barcode number, yaxis is cumulative total # reads
    ax.set_title(f"{type_of_reads} Reads")
    ax.set_xlim((1,max(barsize)))
    ax.set_xscale("log", base=10)
    ax.set_xlabel("Bargroups Ranked by Size") # number of reads per barcode
    ax.set_ylabel("Cumulative Reads Used")
    ax.vlines(vline, 0, 1)
    return

def get_num_16s_reads_for_barcodes(unfiltered_barcode_summary_tsv: str):
    # only 16s reads
    df = pd.read_csv(unfiltered_barcode_summary_tsv, sep="\t", index_col = "Barcode")
    return df["Total # of 16s reads"]

def get_num_reads_for_classified_barcodes(unfiltered_barcode_summary_tsv: str, primers_file: str):
    # 16s and arg reads, but not unclassified reads
    cols_to_use = ["Total # of 16s reads"] + get_arg_names(primers_file)
    df = pd.read_csv(unfiltered_barcode_summary_tsv, sep="\t", index_col = "Barcode")
    # remove unwanted columns
    for col in df.columns:
        if col not in cols_to_use:
            df.drop(col, axis=1, inplace = True)
    return df.sum(axis = 1)

def get_num_reads_for_all_barcodes(b_with_ids: str):
    # all reads: 16s, arg, & unclassified
    num_reads_for_barcodes = []
    with open(b_with_ids, 'r') as f:
        for line in f:
            _16s_reads, arg_reads, unclassified_reads = line.strip().split(": ")[1].split(" |")
            all_reads = _16s_reads.split(", ") + arg_reads.split(", ") + unclassified_reads.split(", ")
            
            # remove empty elements (these are not reads)
            for _ in range(all_reads.count("")):
                all_reads.remove("")
            for _ in range(all_reads.count(" ")):
                all_reads.remove(" ")
            for _ in range(all_reads.count("  ")):
                all_reads.remove("  ")
           
            num_reads_for_barcodes.append(len(all_reads))

    return pd.Series(num_reads_for_barcodes)

# ===============================================================================================

# Helper functions for Primer Balance Figure

def get_ARG_to_16s_ratios(filtered_counts_summary_tsv: str, primers_file: str):
    """
    Calculate the ratio of # of reads of ARG/(ARG + 16s), for each cell with a given ARG
    These ratios will be used in evaluating the primer balance for each ARG - ideally want ratios close to 50/50 
    """

    df_asv_summary = pd.read_csv(filtered_counts_summary_tsv, sep="\t", index_col = "Barcode") # load ASV barcode summary
    arg_names = get_arg_names(primers_file)
    arg_ratios = {}

    for arg in arg_names:
        df_one_arg = df_asv_summary[df_asv_summary[arg] > 0] # get all barcodes that have the given arg 
        if len(df_one_arg) == 0:
            continue # continue if no cell has this ARG

        arg_ratios[arg] = (df_one_arg[arg] / (df_one_arg[arg] + df_one_arg["Total # of 16s reads"])).to_list()

    return arg_ratios


def main():
    parser = argparse.ArgumentParser()
  
    # take input parameters
    parser.add_argument("--use_asvs_str", type=str, required=True)
    parser.add_argument("--unfiltered_barcode_summary_tsv", type=str, required=True)
    parser.add_argument("--final_asv_barcode_summary_tsv", type=str, required=True)
    parser.add_argument("--asv_barcode_summary_no_sub_args_tsv", type=str, required=True)
    parser.add_argument("--primers_file", type=str, required=True)
    parser.add_argument("--b_with_ids", type=str, required=True)
    # 2026-08-10: Route plotting summaries to tmp and final images to figures.
    # Reason: the aggregated taxa-by-target table is a plotting intermediate, not the cell-level result.
    parser.add_argument("--asv_arg_table_tsv", type=str, default="tmp/taxa_target_summary.tsv")
    parser.add_argument("--global_asv_tsv", type=str, required=True)
    parser.add_argument("--asv_arg_figure", type=str, default="figures/taxa_target_table.png")
    parser.add_argument("--barcode_group_size_figure", type=str, default="figures/barcode_group_size_qc.png")
    parser.add_argument("--primer_balance_figure", type=str, default="figures/primer_balance_qc.png")
    parser.add_argument("--first_gene_column_num", type=int, required=True)
    parser.add_argument("--global_mle_tax_tsv", type=str, default="tmp/global_mle_tax.tsv")
    parser.add_argument("--arg_threshold", type=float, default=0.01)
    # 2026-08-10: Keep the CLI default aligned with the 30-cell ASV visualization threshold.
    # Reason: command-line and direct function runs must use the same moderate-sample cutoff.
    parser.add_argument("--min_cells_per_asv", type=int, default=30)
    parser.add_argument("--figure_dpi", type=int, default=300) 
    
    args = parser.parse_args()

    # 2026-08-10: Materialize figure, report, and temporary directories before plotting.
    # Reason: visualization products should not clutter the result root.
    ensure_output_directories(
        args.asv_arg_table_tsv, args.asv_arg_figure,
        args.barcode_group_size_figure, args.primer_balance_figure,
        args.global_mle_tax_tsv)
    
    # make sure input file paths exist
    if not os.path.exists(args.unfiltered_barcode_summary_tsv):
        print(f"❌ Error: input file not found: {args.unfiltered_barcode_summary_tsv}")
        return
    if not os.path.exists(args.final_asv_barcode_summary_tsv):
        print(f"❌ Error: input file not found: {args.final_asv_barcode_summary_tsv}")
        return
    if not os.path.exists(args.asv_barcode_summary_no_sub_args_tsv):
        print(f"❌ Error: input file not found: {args.asv_barcode_summary_no_sub_args_tsv}")
        return
    if not os.path.exists(args.primers_file):
        print(f"❌ Error: input file not found: {args.primers_file}")
        return
    if not os.path.exists(args.b_with_ids):
        print(f"❌ Error: input file not found: {args.b_with_ids}")
        return
    # global_asv_tsv file only needs to exist if we're using ASVs
    if args.use_asvs_str == "yes" and not os.path.exists(args.global_asv_tsv):
        print(f"❌ Error: input file not found: {args.global_asv_tsv}")
        return

    make_figures(
        args.use_asvs_str, args.unfiltered_barcode_summary_tsv, 
        args.final_asv_barcode_summary_tsv, 
        args.asv_barcode_summary_no_sub_args_tsv, args.primers_file, 
        args.b_with_ids, args.asv_arg_table_tsv,  
        args.asv_arg_figure, args.barcode_group_size_figure, args.primer_balance_figure,
        args.first_gene_column_num, args.global_mle_tax_tsv, args.global_asv_tsv,
        args.arg_threshold, args.min_cells_per_asv, args.figure_dpi)


if __name__ == "__main__":
    main()
