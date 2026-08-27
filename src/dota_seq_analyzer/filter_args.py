#!/usr/bin/env python3
import math
import pandas as pd
from helper_functions import get_arg_names, ensure_output_directories
import argparse
import os

def find_poisson_survival_probability(k, lambda_) -> float:
    """
    Determine and return the probability of observing >= k reads,
     given the value of lambda for this specific cell and specific ARG.
     Note that k = the actual observed ARG read count.
    """
    
    if lambda_ == 0:
        return 0
    # Probability of observing >= k reads purely by random background ambient noise
    # P(X >= k) = 1 - P(X <= k-1)
    # P(X <= k-1) = sum_{m=0}^{k-1} (lambda^m * exp(-lambda)) / m!
    sum = 0
    term = math.e**(-1*lambda_); # Initial term where m = 0
    sum += term
    
    m = 1
    while m <= (k - 1):
        term = term * lambda_ / m
        sum += term
        m += 1
    
    p_noise = 1 - sum
    return p_noise # return the probability: P(X >= k | lambda)

def filter_args(
    asv_barcode_summary_tsv: str, primers_file,
    filtered_counts_summary_tsv: str, filtered_binary_summary_tsv: str, stats_filtering_summary_tsv: str, 
    alpha: float = 0.05):
    """
    Remove background ARG noise - that is, cells with too few reads of a given ARG should not be considered to have that ARG.
    Hence nullify these low ARG read count values in the barcode summary.
    Use an algorithm based on the Poisson distribution to determine whether or not a given read count is sufficient,
     based on both the specific cell and the specific ARG (this entails a bidirectional normalization).
    # 2026-08-27: Fix the function docstring terminator.
    # Reason: the extra quote caused a SyntaxError before the filtering stage could run.
    """

    # Note that alpha is the family-wise error rate / significance threshold
    # This controls false positives per cell
        
    # preliminary step
    arg_names = get_arg_names(primers_file)

    # input and output dataframes
    df_original_asv_summary = pd.read_csv(asv_barcode_summary_tsv, sep="\t", index_col = "Barcode")
    df_filtered_counts = pd.DataFrame(df_original_asv_summary)
    df_filtered_binary = pd.DataFrame(df_original_asv_summary)

    df_filtered_counts[arg_names] = df_filtered_counts[arg_names].astype("Int64")
    df_filtered_binary[arg_names] = df_filtered_binary[arg_names].astype("Int64")

    # to be used in statistical summary:
    df_stats = pd.DataFrame(index = arg_names, columns = ["Original", "Filtered Out", "Remaining"])
    df_stats[df_stats.columns] = 0

    # calculate global sequencing scale and aggregate ARG totals
    arg_global_sum = df_original_asv_summary[arg_names].sum(axis=0).to_dict()
    total_16s_global = df_original_asv_summary["Total # of 16s reads"].sum(axis=0)

    # core algorithm; filtering and statistics gathering

    # iterate through each barcode/cell
    for idx in df_original_asv_summary.index.to_list():
        num_16s_reads = df_original_asv_summary["Total # of 16s reads"].loc[idx]

        # iterate through each ARG
        for arg in arg_names:
            k = df_original_asv_summary.loc[idx, arg] # raw number of reads of the specific arg in this cell

            if k > 0:
                # record that this gene was originally positive for this cell
                df_stats.loc[arg, "Original"] += 1

            if (k == 0 or arg_global_sum[arg] == 0):
                # if no reads at all for the specific arg in this cell
                df_filtered_counts.loc[idx, arg] = 0
                df_filtered_binary.loc[idx, arg] = 0

            else:
                # if there are reads, check if # of reads meets min threshold, using Poisson probability
                # Bidirectional normalization to calculate lambda
                r_j = arg_global_sum[arg] / total_16s_global
                lambda_ij = num_16s_reads * r_j
                
                # Execute statistical evaluation to determine probability of observing >=k reads, given the value of lambda
                p_value = find_poisson_survival_probability(k, lambda_ij)

                # Statistical Decision Engine
                if p_value < alpha:
                    df_filtered_counts.loc[idx, arg] = k  # Retain raw quantitative depth
                    df_filtered_binary.loc[idx, arg] = 1  # Convert to positive assignation (1)
                    df_stats.loc[arg, "Remaining"] += 1    # Cell passed the filter
                else:
                    df_filtered_counts.loc[idx, arg] = 0  # Technical noise, clear to 0
                    df_filtered_binary.loc[idx, arg] = 0
                    df_stats.loc[arg, "Filtered Out"] += 1      # Cell was filtered out

    # prepare output stats file
    df_stats["% Retention"] = 0.0
    for arg in df_stats.index.to_list():
        if df_stats.loc[arg, "Original"] == 0:
            df_stats.loc[arg, "% Retention"] = 0.0
        else:
            df_stats.loc[arg, "% Retention"] = df_stats.loc[arg, "Remaining"] / df_stats.loc[arg, "Original"] * 100

    # write to output files
    df_filtered_counts.to_csv(filtered_counts_summary_tsv, sep = "\t", index_label = "Barcode")
    df_filtered_binary.to_csv(filtered_binary_summary_tsv, sep = "\t", index_label = "Barcode")
    df_stats.to_csv(stats_filtering_summary_tsv, sep = "\t", index_label = "ARG")


def main():
    parser = argparse.ArgumentParser()

    # take input parameters
    parser.add_argument("--input_arg_barcode_summary_tsv", type=str, required=True)
    parser.add_argument("--primers_file", type=str, required=True)
    # 2026-08-10: Route ARG-filtering tables to tmp by default.
    # Reason: they are intermediate inputs to subtyping and JSONL export.
    parser.add_argument("--filtered_counts_summary_arg_tsv", type=str, default="tmp/filtered_counts_summary_arg.tsv")
    parser.add_argument("--filtered_binary_summary_arg_tsv", type=str, default="tmp/filtered_binary_summary_arg.tsv")
    parser.add_argument("--stats_filtering_summary_arg_tsv", type=str, default="tmp/stats_filtering_summary_arg.tsv")
    parser.add_argument("--alpha", type=float, default=0.05)

    args = parser.parse_args()

    # 2026-08-10: Materialize the temporary output directory before writing ARG tables.
    # Reason: the default tmp paths must work in a new result directory.
    ensure_output_directories(
        args.filtered_counts_summary_arg_tsv, args.filtered_binary_summary_arg_tsv,
        args.stats_filtering_summary_arg_tsv)
    
    # make sure input file paths exist
    if not os.path.exists(args.input_arg_barcode_summary_tsv):
        print(f"❌ Error: input file not found: {args.input_arg_barcode_summary_tsv}")
        return
    if not os.path.exists(args.primers_file):
        print(f"❌ Error: input file not found: {args.primers_file}")
        return

    filter_args(
        args.input_arg_barcode_summary_tsv, args.primers_file,
        args.filtered_counts_summary_arg_tsv, args.filtered_binary_summary_arg_tsv, 
        args.stats_filtering_summary_arg_tsv, args.alpha)


if __name__ == "__main__":
    main()
