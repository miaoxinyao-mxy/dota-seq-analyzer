#!/usr/bin/env python3
import pandas as pd
import math
from helper_functions import open_maybe_gzip, get_arg_names, get_target_modes, ensure_output_directories
from filter_sub_args import run_sub_arg_denoising_pipeline
from collections import Counter, defaultdict
from typing import List, Dict, Tuple
import os
import argparse
import sys
import json

# 1-based extraction coordinates to accurately truncate and merge R1/R2 reads 
# while safely skipping the 20bp cell barcode sequence at the start of R2, as well as the ~20bp overlap sequence.
# 2026-08-28: Match the sub-locus core endpoints to the ASV core endpoints.
# Reason: both analyses should use the same paired-read sequence region.
R1_START, R1_END = 30, 120
R2_START, R2_END = 70, 120

def create_sub_arg_barcode_summary(
    filtered_counts_summary_arg_tsv: str, 
    b_with_ids: str, arg_packets: str, 
    fwd_fastq: str, rev_fastq: str,
    sub_arg_barcode_summary_tsv: str, primers_file: str, 
    sub_arg_seqs_list: str, 
    filtered_sub_arg_barcode_summary_tsv: str, 
    filtered_stats_cells_per_sub_arg_tsv: str, extra_mle_info_sub_arg_tsv,
    p_match, p_none, p_error, alpha_prior, beta_prior,
    min_confidence, min_noise_reads, noise_cutoff_ratio, include_all_targets=False,
    alpha=0.05, max_shift_sub_arg=2, max_mm_sub_arg=0):

    """
    Main function for creating the sub-ARG barcode summary.
    Note that "sub-ARG" refers to a subtype of the parent ARG.
    This new version of the barcode summary will have ARG read counts
     replaced with either the name of the consensus sub-ARG for that ARG in that cell,
     or the parent ARG name (if unable to confidently classify to the sub-ARG level), 
     or "0" if ARG is not present.
    Hence this sub-ARG barcode summary provides more specific information regarding
     which sub-ARGs are present in each cell.
    """
        
    # note that input barcode summary should have already undergone first round of filtering 
    # (i.e. min 16s reads, contamination, unclassified taxonomy, min barcodes),
    # as well as ARG background noise removal

    #====================================================================================================
    
    # Preliminary steps
    # Read in genes_eligible_for_sub_args
    genes_eligible_for_sub_args = get_genes_eligible_for_sub_args(
        primers_file, include_all_targets)
    
    # 2026-08-10: Keep the sub-ARG gene vector entirely determined by the primer CSV.
    # Reason: filtering must not inject a control gene that is absent from the selected primer panel.

    # Obtaining useful parameter values

    arg_names = get_arg_names(primers_file)
    # Make a list containing the gene numbers of all genes eligible for sub-ARGs
    nums_genes_eligible_for_sub_args = convert_gene_names_to_nums(genes_eligible_for_sub_args, arg_names)

    # Load in the input barcode summary (which still has column names as ARGs)
    # Note this summary has already undergone ARG filtering
    df_arg = pd.read_csv(filtered_counts_summary_arg_tsv, sep="\t", index_col = "Barcode")

    #====================================================================================================

    # Determine all unique sub-ARG sequences for each ARG (post-clustering)
    
    # 1) Find the IDs of all ARG reads corresponding to all filtered barcodes.
    #    Also match each ID to its corresponding gene, to be used in further processing.
    print("A. 1) Matching IDs to genes...")
    ids_to_genes = match_ids_to_genes_bc(b_with_ids, arg_packets, filtered_counts_summary_arg_tsv)

    # 2) Determine all unique sub-ARG sequences for each ARG
    print("A. 2) Determining all unique sub-ARG sequences for each ARG...")
    sub_arg_seqs = get_all_sub_arg_seqs(ids_to_genes, arg_names, fwd_fastq, rev_fastq, nums_genes_eligible_for_sub_args)

    # 3) Cluster the sub-ARG sequences for each ARG.
    #    Account for max shift of 2 between clustered sub-ARG sequences, but no mismatches.
    print("A. 3) Clustering sub-ARG sequences for each ARG...")
    print()
    clustered_sub_arg_seqs = cluster_sub_arg_seqs(arg_names, sub_arg_seqs, nums_genes_eligible_for_sub_args, max_mm_sub_arg, max_shift_sub_arg)

    #====================================================================================================
    
    # Convert clustered_sub_arg_seqs into a more usable form for writing the sub-ARG barcode summary

    # 1) Merge the dicts in clustered_sub_arg_seqs into one dictionary
    print("B. 1) Merging clustered_sub_arg_seqs dicts into one dict...")
    sub_arg_seqs_1_dict = merge_into_one_dict(clustered_sub_arg_seqs)

    # 2) Create a sub-ARG dict sorted by barcode, instead of sub-ARG sequence
    print("B. 2) Sorting sub-ARG dict by barcode, instead of sub-ARG sequence...")
    print()
    sub_arg_seqs_by_barcode = sort_sub_args_by_barcode(df_arg, sub_arg_seqs_1_dict)

    #====================================================================================================

    # Create the sub-ARG barcode summary

    # 1) Write new dataframe barcode summary with sub-ARGs. 
    #    Also write file with the names, cell counts, and core sequences of each sub-ARG.
    print("C. 1) Creating the sub-ARG barcode summary...")
    print()
    df_sub_arg = write_sub_arg_barcode_summary(
        df_arg, genes_eligible_for_sub_args, sub_arg_seqs_by_barcode, arg_names, 
        sub_arg_seqs_list, sub_arg_barcode_summary_tsv, extra_mle_info_sub_arg_tsv,
        p_match, p_none, p_error, alpha_prior, beta_prior,
        min_confidence, min_noise_reads, noise_cutoff_ratio)

    # 2) Filter out sub-ARGs with too few cells associated with them, and obtain final list of sub-ARGs to be used
    print("C. 2) Filtering sub-ARGs, and obtaining final names list of sub-ARGs to be used...")
    print()
    final_sub_arg_names = run_sub_arg_denoising_pipeline(
        sub_arg_seqs_list, alpha, filtered_stats_cells_per_sub_arg_tsv)

    # 3) Replace filtered-out sub-ARGs with their respective parent ARG.
    print("C. 3) Modifying barcode summary based on filtered sub-ARGs...")
    print()
    modify_sub_arg_barcode_summary(
        df_sub_arg, final_sub_arg_names, filtered_sub_arg_barcode_summary_tsv,
        genes_eligible_for_sub_args)


def write_sub_arg_barcode_summary(
    df_arg, genes_eligible_for_sub_args, sub_arg_seqs_by_barcode, arg_names, 
    sub_arg_seqs_list, sub_arg_barcode_summary_tsv, extra_mle_info_sub_arg_tsv,
    p_match, p_none, p_error, alpha_prior, beta_prior,
    min_confidence, min_noise_reads, noise_cutoff_ratio):
    """
    Update the barcode summary to include sub-ARG information, and write to a TSV.
    Specifically, replace each numeric ARG read count with either the consensus sub-ARG name,
     parent ARG name (if either a consensus sub-ARG could not be confidently identified, or if the ARG is not eligible for sub-ARGs),
     or "0" if ARG is not present in that cell.
    Arbitrarily name the observed sub-ARGs in order of frequency, and record the sub-ARG sequence for each arbitrary name in a text file.
    In this names-to-sequences file, also include the cell count for each sub-ARG -> this will be used when filtering out sub-ARGs later on.
    When using the MLE consensus algorithm to determine consensus sub-ARGs, record supplementary info (e.g. contamination level) in a separate TSV.
    Return the updated sub-ARG barcode summary as a pandas dataframe.
    """

    # copy the original ARG-based barcode summary to the new sub-ARG-based barcode summary dataframe
    df_sub_arg = pd.DataFrame(df_arg)
    df_sub_arg[arg_names] = df_sub_arg[arg_names].astype(str)

    # make a new dataframe to record supplementary/extra information
    # this info will be generated when determining the consensus sub-ARG for each relevant ARG in each cell, using the MLE consensus algorithm
    barcodes = df_arg.index
    mle_output_cols = ["Predicted_sub-ARG", "Confidence", "Contamination", "Total_#_of_ARG_reads", "Technical_noise_count"]
    cols_for_mle_df = []
    for gene in genes_eligible_for_sub_args:
        for mle_col_name in mle_output_cols:
            cols_for_mle_df.append(f"{gene}: {mle_col_name}")
    df_extra_mle_info = pd.DataFrame(index = barcodes, columns = cols_for_mle_df)
    df_extra_mle_info[:] = "-"

    # iterate through each barcode/cell
    for bc in barcodes:
        # for this specific barcode, obtain all distinct core sequences for all ARGs
        # each sequence will also have the number of the ARG it corresponds to,
        # and the # of reads having that sequence
        read_counts_all_args = sub_arg_seqs_by_barcode[bc]

        # iterate through each ARG that is eligible to have sub-ARGs
        for arg in genes_eligible_for_sub_args:

            gene_num = arg_names.index(arg) # get the number corresponding to the current ARG

            # continue to next iteration, if cell does not contain this ARG
            if df_arg.loc[bc, arg] == 0:
                continue

            # collect all distinct core sequences for the given ARG in the given cell, 
            # along with the count of # of reads having each sequence,
            # into the read_counts dictionary
            read_counts = {}
            for seq, count in read_counts_all_args.items():
                # note that seq is in the form of e.g. "ATCGATCG_2", 
                # where "ATCGATCG" is the core sequence, and 2 refers to the ARG number
                if int(seq.split("_")[1]) == gene_num:
                    read_counts[seq.split("_")[0]] = count

            # determine (if possible) one consensus sub-ARG for this ARG in this cell, using the MLE consensus algorithm
            final_sub_arg, total_reads, technical_noise_count, last_valid_confidence, bayesian_contamination_mean \
            = sub_arg_parse_and_analyze_perfect_corrected(
            read_counts, p_match, p_none, p_error, alpha_prior, beta_prior,
            min_confidence, min_noise_reads, noise_cutoff_ratio)

            # record info in sub-ARG barcode summary, based on whether or not a consensus sub-ARG was determined
            if final_sub_arg is not None:
                # if a consensus sub-ARG was determined
                # sub-ARG barcode summary: replace ARG read count with the sequence of the consensus sub-ARG
                df_sub_arg.loc[bc, arg] = final_sub_arg
                # extra MLE info table: record info outputted by MLE consensus algorithm
                df_extra_mle_info.loc[bc, f"{arg}: Predicted_sub-ARG"] = final_sub_arg
                df_extra_mle_info.loc[bc, f"{arg}: Confidence"] = last_valid_confidence
                df_extra_mle_info.loc[bc, f"{arg}: Contamination"] = bayesian_contamination_mean
                df_extra_mle_info.loc[bc, f"{arg}: Total_#_of_ARG_reads"] = total_reads
                df_extra_mle_info.loc[bc, f"{arg}: Technical_noise_count"] = technical_noise_count
            else:
                # if a consensus sub-ARG could not be confidently determined
                # sub-ARG barcode summary: replace ARG read count with "[ARG]_parent", where [ARG] is replaced by the name of the current ARG
                df_sub_arg.loc[bc, arg] = f"{arg}_parent"

    # assign standardized names to each consensus sub-locus sequence and document these mappings
    # note this step is done after consensus sub-ARGs have been identified for all cells,
    # so that the naming of these sub-ARG sequences can be organized by frequency of the sequences on a global scale
    with open(sub_arg_seqs_list, 'w') as f:

        f.write("Sub-ARG_Arbitrary_Name\tCell_count\tCore_sequence\n")

        # name each final sub-locus sequence in descending global frequency order
        # replace sequences with their new names in the barcode summary
        # also record these names-to-sequences in a text file
        for arg in genes_eligible_for_sub_args:
            # obtain all consensus sub-ARG sequences for this ARG, across all cells
            final_sub_arg_seqs = df_sub_arg.loc[ \
                (df_sub_arg[arg] != f"{arg}_parent") & \
                (df_sub_arg[arg] != "0") & \
                (df_sub_arg[arg] != "0.0"), \
                arg].to_list()

            # name sub-ARG sequences from most to least common, on a global scale
            # e.g. TEM_seq_1 would have more cells associated with it than TEM_seq_2
            ranked_final_sub_arg_seqs = Counter(final_sub_arg_seqs).most_common()

            for i, (seq, cell_count) in enumerate(ranked_final_sub_arg_seqs, start=1):
                sub_arg_name = f"{arg}_seq_{i}" # 2026-08-28: Use parent_seq_N identifiers for generated sequence clusters.
                # Reason: distinguish pipeline-generated clusters from pre-existing biological variant names.
                f.write(f"{sub_arg_name}\t{cell_count}\t{seq}\n") # record the sequence associated with this new arbitrary sub-ARG name in a text file
                df_sub_arg.loc[df_sub_arg[arg] == seq, arg] = sub_arg_name # update the sub-ARG barcode summary to replace sub-ARG sequences with their new names

                df_extra_mle_info.loc[df_extra_mle_info[f"{arg}: Predicted_sub-ARG"] == seq, \
                    f"{arg}: Predicted_sub-ARG"] = sub_arg_name     # update the extra MLE info table to replace sub-ARG sequences with their new names

    # update sub-ARG barcode summary for ARGs that are not eligible for sub-ARGs
    # specifically, replace numeric ARG read counts with either the ARG name or "0"
    for bc in df_sub_arg.index:
        for arg in arg_names:
            # skip over the sub-ARG columns, because we already handled those
            if arg in genes_eligible_for_sub_args:
                continue
            # replace numeric read counts with either the ARG name (if present), or "0" (if absent)
            try:
                num = int(float(df_sub_arg.loc[bc, arg]))
                if num != 0:
                    df_sub_arg.loc[bc, arg] = arg
            except:
                assert False # original values should be integer read counts, so return assertion error if not

    # write dataframes to output TSV files
    df_sub_arg.to_csv(sub_arg_barcode_summary_tsv, sep = "\t", index_label = "Barcode")
    df_extra_mle_info.to_csv(extra_mle_info_sub_arg_tsv, sep = "\t", index_label = "Barcode")

    return df_sub_arg

def modify_sub_arg_barcode_summary(
    df_sub_arg, final_sub_arg_names, filtered_sub_arg_barcode_summary_tsv,
    genes_eligible_for_sub_args):
    """Replace filtered-out sub-ARG names with their parent ARG names, in the barcode summary"""
        
    # 2026-08-10: Apply one consistent parent-gene fallback to primer-selected genes.
    # Reason: there is no longer a user-visible or forced baseline gene.

    # for cells containing sub-ARGs that were filtered out, replace the sub-ARG name with the parent ARG name
    for arg in genes_eligible_for_sub_args:
        df_sub_arg.loc[~df_sub_arg[arg].isin(["0", f"{arg}_parent"] + final_sub_arg_names), arg] = f"{arg}_parent"

    # write to TSV
    df_sub_arg.to_csv(filtered_sub_arg_barcode_summary_tsv, sep = "\t", index_label = "Barcode")


def match_ids_to_genes_bc(
    b_with_ids: str, arg_packets: str, 
    filtered_barcode_summary_tsv: str) -> Dict[str, int]:
    """
    Goal: obtain information for all ARG read IDs related to the current filtered barcode summary.
    Create a dictionary of the form {"ID": "gene_num|barcode", ...}

    This dictionary will contain all IDs that are both 1) ARG reads, 
     and 2) ID corresponds to one of the barcodes in the input filtered_barcode_summary_tsv 
     (i.e. some barcode filtering has already been done - this ensures we avoid ARG read with sequences that are too short,
      as that would affect later processing, specifically the extract_core function).

    This dictionary will be used in downstream processing.
    """

    # obtain list of all current barcodes (i.e. excludes filtered-out barcodes)
    filtered_bcs = pd.read_csv(filtered_barcode_summary_tsv, sep="\t", index_col = "Barcode").index

    # compile list of all ARG reads' IDs with their corresponding barcodes, with each element being a string in the form of "ID|Barcode"
    filtered_arg_ids_with_bs = []
    with open(b_with_ids, 'r') as b_with_ids_file:
        for line in b_with_ids_file:
            bc, ids = line.strip().split(": ")
            if bc in filtered_bcs: # barcode must be in current filtered barcode list (meaning not already filtered out)
                arg_ids_field = ids.split(" | ")[1]
                # 2026-08-10: Detect ARG read IDs from field content instead of an instrument-specific prefix.
                # Reason: valid FASTQ IDs do not necessarily contain "SH0", and must not be silently discarded.
                arg_ids = [
                    read_id.strip(" |")
                    for read_id in arg_ids_field.split(",")
                    if read_id.strip(" |")
                ]
                arg_ids_with_bs = [f"{read_id}|{bc}" for read_id in arg_ids]
                filtered_arg_ids_with_bs.extend(arg_ids_with_bs)

    # extract list of the IDs of all ARG reads corresponding to barcodes in the current filtered barcode list
    filtered_arg_ids = [id_with_b.split("|")[0] for id_with_b in filtered_arg_ids_with_bs]

    # compile a dictionary mapping each ARG read's ID to its corresponding gene number and barcode
    # this ids_to_genes_bc dictionary will be in the form of {"ID": "gene_num|barcode", ...}
    ids_to_genes_bc = {}
    with open(arg_packets, 'r') as arg_packet_file:
        # iterate through each ARG read packet
        for line in arg_packet_file:
            packet = json.loads(line.strip())
            id = packet["ID"] # get ID
            if id in filtered_arg_ids: # only use IDs that are in our filtered IDs list
                bc = filtered_arg_ids_with_bs[filtered_arg_ids.index(id)].split("|")[1] # get barcode
                gene_num = packet["gene"].index(1) # get gene number
                ids_to_genes_bc[id] = f"{gene_num}|{bc}" # add to the dictionary
    return ids_to_genes_bc

def get_all_sub_arg_seqs(
    ids_to_genes: Dict, arg_names: List[str], fwd_fastq: str, rev_fastq: str, 
    nums_genes_eligible_for_sub_args: List[int]) -> List[Dict]:
    """
    Obtain all sub-ARG sequences present, for all ARGs that are eligible to have sub-ARGs.

    Returned value is a list of dictionaries, where each dictionary corresponds to one ARG.
    Each dictionary is of the form: (e.g. of TEM as being the parent ARG for this dictionary)
    {"Sub-ARG seq #1": ["ID#5|Barcode#30", "ID#8|Barcode#7"],   # corresponds to e.g. TEM_seq_1
     "Sub-ARG seq #2": ["ID#10|Barcode#1"], ...}                # corresponds to e.g. TEM_seq_2
    Note that all values chosen for this example are random - the goal above was just to show
    the sub_arg_seqs formatting.
    Each ID corresponding to a specific sub-ARG sequence are added to that sub-ARG's list
    (e.g. ID#5 and ID#8 would correspond to TEM_seq_1 in the example above.
    The ID names are from the fastq file. Barcodes are included with the IDs in the format of
    "ID|Barcode", to faciliate later downstream processing.
    Note also that "Sub-ARG seq #1" would be in the format outputted by the extract_core() function;
    that is, "ATCG...|TACC..." (with the specific sequences provided here being just an example).

    Return a list of dictionaries mapping each sub-ARG sequence to all its relevant IDs,
    while tracking parent ARGs (for the sub-ARGs) and barcodes (for the IDs).
    """

    sub_arg_seqs = [{} for _ in range(len(arg_names))] 

    # 2026-09-02: Count eligible ARG reads that cannot provide a complete core.
    # Reason: extract_core() returns None for an R1 or R2 shorter than the required 120 bp.
    skipped_short_reads = 0

    with open_maybe_gzip(fwd_fastq, 'r') as fwd, open_maybe_gzip(rev_fastq, 'r') as rev:
        # read in the read_ID line
        f_line = fwd.readline()
        r_line = rev.readline()

        # iterate through each read line in the fastq files
        i = 0
        while (f_line != "") and (r_line != ""):

            # parsing
            
            id_f = f_line.strip().split(" ")[0].strip("@")
            id_r = r_line.strip().split(" ")[0].strip("@")
            # check that IDs for this line are equal - necessary for this analysis
            assert id_f == id_r, f"ID from forward and reverse fastq files do not match on line {i*4 + 1}"
            # sequence parsing
            f_seq, r_seq = fwd.readline().strip(), rev.readline().strip()
            # prepare for reading in the next read_ID line
            for _ in range(3):
                f_line = fwd.readline()
                r_line = rev.readline()

            # if ID is in the filtered ARG IDs list (and thus in the ids_to_genes dictionary),
            # then obtain its gene and barcode info
            gene_num_and_bc = ids_to_genes.get(id_f)
            if gene_num_and_bc != None:
                # parsing
                gene_num, bc = gene_num_and_bc.split("|")
                gene_num = int(gene_num)

                # skip over ARGs that are not eligible to have sub-ARGs
                if gene_num not in nums_genes_eligible_for_sub_args:
                    i += 1 # ensure counter still increments, before "continue" keyword
                    continue

                # obtain core of the R1 & R2 sequences for this read
                # then add this to the appropriate position in the sub_arg_seqs dictionary
                seq_pair = extract_core(f_seq, r_seq)
                if seq_pair is None:
                    # 2026-09-02: Do not pass incomplete cores into sequence clustering.
                    # Reason: clustering expects an R1|R2 string and cannot process None.
                    skipped_short_reads += 1
                    i += 1
                    continue
                if seq_pair in sub_arg_seqs[gene_num].keys():
                    sub_arg_seqs[gene_num][seq_pair].append(f"{id_f}|{bc}")
                else:
                    sub_arg_seqs[gene_num][seq_pair] = [f"{id_f}|{bc}"]
            i += 1

    print(f"Skipped {skipped_short_reads} eligible ARG reads with R1 or R2 shorter than {R1_END} bp.")

    return sub_arg_seqs

def cluster_sub_arg_seqs(
    arg_names: List[str], sub_arg_seqs: List[Dict], nums_genes_eligible_for_sub_args: List[int], 
    max_mm_sub_arg: int, max_shift_sub_arg: int) -> List[Dict]:
    """
    Cluster the sub-ARG sequences, to create a whitelist of sub-ARG sequences for each ARG.
    Sub-ARG sequences are clustered based on the given max mismatch (should be zero) and max shift,
     along with the fact that the cluster boundary should be maintained.
    Returns the modified sub_arg_seqs list (now called clustered_sub_arg_seqs).
    """

    clustered_sub_arg_seqs = [{} for _ in range(len(arg_names))]

    # iterate through each ARG
    for gene_num in range(len(sub_arg_seqs)):

        # skip over ARGs that are not eligible to have sub-ARGs
        if gene_num not in nums_genes_eligible_for_sub_args:
            continue

        # count the number of reads associated with each distinct sub-ARG sequence
        # then sort sub-ARG sequences from most to least common
        counts = {}
        for seq, ids in sub_arg_seqs[gene_num].items():
            counts[seq] = len(ids)
        sorted_sub_seqs = [k for k, v in Counter(counts).most_common()]

        # Iterate through sub-ARG sequences, from most to least common.
        # Reason: Dominant sequences are identified as we go through this list - hence more dominant sequences will be identified near the start of the list.
        #         Note that the dominant sequences are those representing their respective sub-ARG clusters.
        #         We want the dominant sequences to be the most common sequence in their respective clusters - 
        #         combined with the fact that more dominant sequences will be identified near the list's beginning, 
        #         this is why we must sort the sub-ARG sequence list by size (i.e. # of associated reads) before iterating through it.
        for seq_s in sorted_sub_seqs:
            in_dominant_seqs = False
            # iterate through all dominant sub-ARG sequences identified so far for that ARG
            for seq_d in clustered_sub_arg_seqs[gene_num]: 
                d = paired_semi_global_distance(seq_s, seq_d, max_shift=max_shift_sub_arg) # calculate distance between current sequence and given dominant sequence
                # if current sequence can be clustered with one of the dominant sequences, then merge their read ID lists
                if d <= max_mm_sub_arg and within_cluster_boundary(clustered_sub_arg_seqs[gene_num], seq_s, max_mm_sub_arg, max_shift_sub_arg):
                    in_dominant_seqs = True
                    clustered_sub_arg_seqs[gene_num][seq_d].extend(sub_arg_seqs[gene_num][seq_s]) # merge the read ID lists for that cluster
                    seq_added = True
                    break # similar sequence already in dominant sequences list
                    
            # if current sequence cannot be clustered with any dominant sequence, 
            # then add it to clustered_sub_arg_seqs as a new dominant sequence
            if not in_dominant_seqs:
                clustered_sub_arg_seqs[gene_num][seq_s] = sub_arg_seqs[gene_num][seq_s]
                seq_added = True

            # sequence should have been added to clustered_sub_arg_seqs already; if not, return error
            assert seq_added == True, f"Sequence {seq_s} not matched to any sequence in the dominant sequences whitelist for ARG # {gene_num}"

    return clustered_sub_arg_seqs

def extract_core(
    r1: str, r2: str) -> str:
    """Slices and concatenates the predefined hypervariable core regions from R1 and R2."""
    if len(r1) < R1_END or len(r2) < R2_END:
        return None
    return r1[R1_START:R1_END] + "|" + r2[R2_START:R2_END]

def merge_into_one_dict(
    clustered_sub_arg_seqs: List[Dict]) -> Dict:
    """
    Merge all sub-ARG sequences (across all ARGs) into a single dictionary, to make it easier to organize info 
     on a per-barcode basis later on, thus faciliating downstream per-barcode analysis.
    To ensure we can still identify the parent ARG for each sub-ARG sequence, append the gene number to the end of each sequence. 
    Return this single merged dictionary.
    """

    # preliminary steps
    sub_arg_seqs_1_dict = {}
    total_num_sub_arg_seqs = 0

    # iterate through each ARG (note that i is the gene number)
    for i in range(len(clustered_sub_arg_seqs)):
        new_one_arg_dict = {}
        # add each sequence with its corresponding IDs to the new dictionary
        # also append the gene number to the end of each sequence, so that the gene can be identified, 
        # even when all sub-ARG sequences (across all ARGs) are compiled into the same dictionary
        for seq, ids in clustered_sub_arg_seqs[i].items():
            new_one_arg_dict[f"{seq}_{i}"] = ids

        # add all elements from the temporary dictionary for this gene into the single final dictionary
        sub_arg_seqs_1_dict.update(new_one_arg_dict)
        total_num_sub_arg_seqs += len(new_one_arg_dict)        

    # check that no sub-ARG sequences were lost when updating dict
    assert len(sub_arg_seqs_1_dict) == total_num_sub_arg_seqs 

    return sub_arg_seqs_1_dict

def sort_sub_args_by_barcode(
    df_arg, sub_arg_seqs_1_dict: Dict) -> Dict:
    """
    Organize sub-ARG sequence information on a per-barcode basis.
    This will make it easier later on to compile the sub-ARG barcode summary (which uses a per-barcode basis as well).
    Return a dictionary containing this sub-ARG sequence info, sorted by barcode.
    """
    # set up a dictionary, with the keys being all the current filtered barcodes
    filtered_bcs = df_arg.index
    sub_arg_seqs_by_barcode = dict.fromkeys(filtered_bcs)
    for bc in sub_arg_seqs_by_barcode:
        sub_arg_seqs_by_barcode[bc] = {}

    # organize info on a per-barcode basis
    # the sub_arg_seqs_by_barcode dictionary should be in the form of 
    # {barcode1: {seq1: 5, seq2: 1, ...}, 
    #  barcode2: {seq3: 4, ...}, ...}
    # where the numbers corresponding to each sequence are the # of reads associated with that sequence
    for seq, ids_with_bs in sub_arg_seqs_1_dict.items():
        for id_and_b in ids_with_bs:
            bc = id_and_b.split("|")[1]
            if seq in sub_arg_seqs_by_barcode[bc]:
                sub_arg_seqs_by_barcode[bc][seq] += 1
            else:
                sub_arg_seqs_by_barcode[bc][seq] = 1

    return sub_arg_seqs_by_barcode

def get_genes_eligible_for_sub_args(
    primers_file: str, include_all_targets: bool = False) -> List[str]:
    """
    Extract information from primers file to form a list of all genes eligible to have sub-ARGs.
    The ARGs in the primers file are classified as either "single" or "family".
    Those classified as "family" are eligible to have sub-ARGs.
    Returns a list of eligible genes.
    """
    # 2026-08-11: Reconstruct sequences for PV targets or every target when a reference is supplied.
    # Reason: blank targets need only detection unless reference annotation was explicitly requested.
    target_modes = get_target_modes(primers_file)
    if include_all_targets:
        return list(target_modes)
    return [target for target, mode in target_modes.items() if mode == "ssr"]

def convert_gene_names_to_nums(
    selected_gene_names: List[str], arg_names: List[str]) -> List[int]:
    """
    Convert a list of ARG names to their corresponding zero-indexed numbers.
    To be used for the genes_eligible_for_sub_args list.
    """
    nums_selected_genes = [] 
    for arg in selected_gene_names:
        nums_selected_genes.append(arg_names.index(arg))
    return nums_selected_genes

def semi_global_distance(
    a: str, b: str, max_shift=2) -> int:
    """
    Calculates the minimum mismatch distance between two sequences allowing 
    a sliding window shift to adjust for potential sequencing indels/offsets.
    """

    best = 999999
    for shift in range(-max_shift, max_shift + 1):
        if shift >= 0:
            a_sub = a[shift:]
            b_sub = b[:len(a_sub)]
        else:
            b_sub = b[-shift:]
            a_sub = a[:len(b_sub)]

        n = min(len(a_sub), len(b_sub))
        if n == 0:
            continue

        a_sub = a_sub[:n]
        b_sub = b_sub[:n]

        mismatches = sum(x != y for x, y in zip(a_sub, b_sub))
        best = min(best, mismatches)
    return best

def paired_semi_global_distance(a, b, max_shift=2):
    """Compare paired R1|R2 cores while allowing each read its own shift."""
    a_r1, a_r2 = a.split("|", 1)
    b_r1, b_r2 = b.split("|", 1)
    # 2026-08-28: Compare R1 and R2 independently instead of sharing one shift.
    # Reason: each read can have a valid indexing/length shift independent of its mate.
    return max(
        semi_global_distance(a_r1, b_r1, max_shift=max_shift),
        semi_global_distance(a_r2, b_r2, max_shift=max_shift),
    )

def within_cluster_boundary(
    candidate_seqs: List[str], new_seq: str, 
    max_dist=0, max_shift_sub_arg=2) -> bool:
    """
    Check whether adding new_seq keeps the whole cluster boundary <= max_dist.
    """
    for old_seq in candidate_seqs:
        d = paired_semi_global_distance(new_seq, old_seq, max_shift=max_shift_sub_arg)
        if d > max_dist:
            return False
    return True

def sub_arg_parse_and_analyze_perfect_corrected(
    read_counts: Dict[str, int],
    p_match = 0.90,
    p_none = 0.09,
    p_error = 0.01,
    alpha_prior = 1.0,
    beta_prior = 9.0,
    min_confidence = 0.95,
    min_noise_reads = 2,
    noise_cutoff_ratio = 0.05
):
    """
    Determine whether to classify to the parent ARG or sub-ARG level (and for the latter, determine the consensus sub-ARG).
    Input is the list of all observed sub-ARG sequences with their respective counts, for the given ARG in the given cell.
    Return either the sequence of the consensus sub-ARG; or, if only classified to the ARG level, return None.
    Also return supplementary information (e.g. contamination level).

    The algorithm works in three broad stages:
    1. Use maximum likelihood to choose the dominant sub-ARG sequence.
    2. Only go to sub-ARG rank level if the best sub-ARG assignment has sufficient confidence; 
       otherwise, revert to parent ARG rank level.
    3. Treat small off-path sub-ARGs as technical noise and larger off-path sub-ARGs as
       real contamination before computing a Bayesian contamination estimate.
    """

    total_reads = sum(read_counts.values())

    # ---------------- 1. Stepwise adaptive MLE sub-ARG path inference ----------------
    final_path = ""
    last_valid_confidence = 1.0

    # 🛠️ 【微调点 1】：允许 None 成为合法候选者，不再用 is not None 过滤它
    # 这样当大多数 reads 走到某一层断掉变成 None 时，None 作为一个群体能抱团对抗噪声
    lvl_candidates = read_counts.keys()

    mle_scores = {}
    for cand in lvl_candidates:
        log_likelihood = 0.0
        for r_val, count in read_counts.items():
            
            # 🛠️ 【微调点 2】：针对 None 候选者，设计合理的条件似然度打分
            # 候选是具体菌名（如 Vibrio），维持原来的传统打分逻辑不变
            if r_val == cand:
                for _ in range(count): # account for # of reads with that r_val
                    log_likelihood += math.log(p_match)
            else:
                for _ in range(count):
                    log_likelihood += math.log(p_error)
                    
        mle_scores[cand] = log_likelihood

    # Convert log-likelihoods to relative probabilities with the
    # log-sum-exp trick, avoiding overflow or underflow on many reads.
    max_log = max(mle_scores.values())
    exp_scores = {
        cand: math.exp(score - max_log) for cand, score in mle_scores.items()
    }
    sum_exp = sum(exp_scores.values())

    best_cand = max(mle_scores, key=mle_scores.get)
    confidence = exp_scores[best_cand] / sum_exp

    # If confidence is high enough, accept sub-ARG rank with the best candidate sub-ARG sequence.
    # If not, then accept ARG rank (don't move deeper to sub-ARG rank).
    if confidence >= min_confidence:
        final_path = best_cand
        last_valid_confidence = confidence
    else:
        final_path = None

    # ---------------- 2. Separate technical noise (typing errors) from real contamination ----------------

    # set these values to '-' (equivalent of N/A in the extra MLE info table) by default
    technical_noise_count = "-"
    bayesian_contamination_mean = "-"

    # only calculate contamination level if classified to the sub-ARG rank level
    # reason: if classified to ARG rank level, there should be 100% confidence and no contamination
    if final_path is not None:

        match_reads_count = 0
        real_contamination_count = 0
        technical_noise_count = 0

        for bug, count in read_counts.items():
            if bug == final_path:
                match_reads_count += count
            else:
                bug_ratio = count / total_reads
                if count < min_noise_reads or bug_ratio < noise_cutoff_ratio:
                    technical_noise_count += count
                else:
                    real_contamination_count += count

        # ---------------- 3. Calculate the corrected Bayesian posterior contamination rate ----------------
        corrected_match_count = match_reads_count + technical_noise_count

        alpha_post = alpha_prior + real_contamination_count
        beta_post = beta_prior + corrected_match_count
        bayesian_contamination_mean = alpha_post / (alpha_post + beta_post)

    return final_path, total_reads, technical_noise_count, last_valid_confidence, bayesian_contamination_mean


def main():
    parser = argparse.ArgumentParser()
  
    # take input parameters
    # 2026-08-10: Route sub-ARG working tables to tmp by default.
    # Reason: only the exported JSONL and selected reports are final products.
    parser.add_argument("--filtered_counts_summary_arg_tsv", type=str, default="tmp/filtered_counts_summary_arg.tsv")
    parser.add_argument("--b_with_ids", type=str, required=True)
    parser.add_argument("--arg_packets", type=str, required=True)
    # 2026-08-10: Expose sub-ARG inputs as R1/R2 in the public CLI.
    # Reason: sequencing inputs are named R1 and R2, not forward/reverse.
    parser.add_argument("--r1_fastq", type=str, required=True)
    parser.add_argument("--r2_fastq", type=str, required=True)
    parser.add_argument("--sub_arg_barcode_summary_tsv", type=str, default="tmp/sub_arg_barcode_summary.tsv")
    parser.add_argument("--primers_file", type=str, required=True)
    parser.add_argument("--sub_arg_seqs_list", type=str, default="tmp/sub_arg_seqs_list.txt")
    parser.add_argument("--filtered_sub_arg_barcode_summary_tsv", type=str, default="reports/cell_target_matrix.tsv")
    parser.add_argument("--filtered_stats_cells_per_sub_arg_tsv", type=str, default="tmp/filtered_stats_cells_per_sub_arg.tsv")
    parser.add_argument("--alpha", type=float, default=0.05)
    parser.add_argument("--max_shift_sub_arg", type=int, default=2)
    parser.add_argument("--max_mm_sub_arg", type=int, default=0)
    parser.add_argument("--extra_mle_info_sub_arg_tsv", type=str, default="tmp/extra_mle_info_sub_arg.tsv")
    
    parser.add_argument("--p_match", type=float, default=0.90)
    parser.add_argument("--p_none", type=float, default=0.09)
    parser.add_argument("--p_error", type=float, default=0.01)
    parser.add_argument("--alpha_prior", type=float, default=1.0)
    parser.add_argument("--beta_prior", type=float, default=9.0)
    parser.add_argument("--min_confidence", type=float, default=0.95)
    parser.add_argument("--min_noise_reads", type=int, default=2)
    parser.add_argument("--noise_cutoff_ratio", type=float, default=0.05)
    # 2026-08-10: Permit the CLI to reconstruct every target for optional reference matching.
    # Reason: reference availability is independent of the per-target PV mode.
    parser.add_argument("--include_all_targets", action="store_true")


    args = parser.parse_args()

    # 2026-08-10: Materialize temporary and final report directories before writing sub-ARG tables.
    # Reason: the cell-by-target matrix is a report while working tables remain under tmp.
    ensure_output_directories(
        args.sub_arg_barcode_summary_tsv, args.sub_arg_seqs_list,
        args.filtered_sub_arg_barcode_summary_tsv,
        args.filtered_stats_cells_per_sub_arg_tsv, args.extra_mle_info_sub_arg_tsv)
    
    # make sure input file paths exist
    if not os.path.exists(args.filtered_counts_summary_arg_tsv):
        print(f"❌ Error: input file not found: {args.filtered_counts_summary_arg_tsv}")
        sys.exit(1)
    if not os.path.exists(args.r1_fastq):
        print(f"❌ Error: input file not found: {args.r1_fastq}")
        sys.exit(1)
    if not os.path.exists(args.r2_fastq):
        print(f"❌ Error: input file not found: {args.r2_fastq}")
        sys.exit(1)
    if not os.path.exists(args.primers_file):
        print(f"❌ Error: input file not found: {args.primers_file}")
        sys.exit(1)
    if not os.path.exists(args.b_with_ids):
        print(f"❌ Error: input file not found: {args.b_with_ids}")
        sys.exit(1)
    if not os.path.exists(args.arg_packets):
        print(f"❌ Error: input file not found: {args.arg_packets}")
        sys.exit(1)


    create_sub_arg_barcode_summary(
        args.filtered_counts_summary_arg_tsv, 
        args.b_with_ids, args.arg_packets, 
        args.r1_fastq, args.r2_fastq,
        args.sub_arg_barcode_summary_tsv, args.primers_file, 
        args.sub_arg_seqs_list, 
        args.filtered_sub_arg_barcode_summary_tsv, 
        args.filtered_stats_cells_per_sub_arg_tsv,
        args.extra_mle_info_sub_arg_tsv,
        args.p_match, args.p_none, args.p_error, args.alpha_prior, args.beta_prior,
        args.min_confidence, args.min_noise_reads, args.noise_cutoff_ratio,
        args.include_all_targets,
        args.alpha, args.max_shift_sub_arg, args.max_mm_sub_arg)

if __name__ == "__main__":
    main()
