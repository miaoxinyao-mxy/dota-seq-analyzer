import pandas as pd
import math
from helper_functions import open_maybe_gzip, get_arg_names, ensure_output_directories
from filter_sub_args import run_sub_arg_denoising_pipeline
from collections import Counter, defaultdict
from typing import List, Dict, Tuple
import os
import argparse
import json

# 1-based extraction coordinates to accurately truncate and merge R1/R2 reads 
# while safely skipping the 20bp cell barcode sequence at the start of R2, as well as the ~20bp overlap sequence.
R1_START, R1_END = 30, 110
R2_START, R2_END = 70, 110 

def create_sub_arg_barcode_summary(
    filtered_counts_summary_arg_tsv: str, 
    b_with_ids: str, arg_packets: str, 
    fwd_fastq: str, rev_fastq: str,
    sub_arg_barcode_summary_tsv: str, primers_file: str, 
    sub_arg_seqs_list: str, 
    filtered_sub_arg_barcode_summary_tsv: str, 
    baseline_gene: str, filtered_stats_cells_per_sub_arg_tsv: str, extra_mle_info_sub_arg_tsv,
    p_match, p_none, p_error, alpha_prior, beta_prior,
    min_confidence, min_noise_reads, noise_cutoff_ratio,
    alpha=0.05, max_shift_sub_arg=2, max_mm_sub_arg=0):

    # note that input barcode summary should have already undergone first round of filtering 
    # (i.e. min 16s reads, contamination, unclassified taxonomy, min barcodes),
    # as well as ARG filtering

    #====================================================================================================
    
    # Preliminary steps
    # Read in genes_eligible_for_sub_args
    genes_eligible_for_sub_args = get_genes_eligible_for_sub_args(primers_file)
    
    # Adjusting genes_eligible_for_sub_args to include baseline gene, if not already included

    # Necessary so that baseline gene can be used in filtering sub-ARGs
    # However, if baseline gene is not actually eligible for sub-ARGs, 
    # it will be treated as such in the final barcode summary - so its inclusion in genes_eligible_for_sub_args will be temporary
    if baseline_gene not in genes_eligible_for_sub_args:
        genes_eligible_for_sub_args.append(baseline_gene)
        is_baseline_gene_eligible_for_sub_args = False
    else:
        is_baseline_gene_eligible_for_sub_args = True

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
    final_sub_arg_names = run_sub_arg_denoising_pipeline(sub_arg_seqs_list, 
                        baseline_gene, alpha, filtered_stats_cells_per_sub_arg_tsv)

    # 3) Modify barcode summary - specifically, account for whether or not the baseline gene is actually
    #    eligible for sub-ARGs, and replace any filtered-out sub-ARGs with their respective parent ARG
    print("C. 3) Modifying barcode summary based on filtered sub-ARGs...")
    print()
    modify_sub_arg_barcode_summary(
        df_sub_arg, final_sub_arg_names, filtered_sub_arg_barcode_summary_tsv,
        baseline_gene, is_baseline_gene_eligible_for_sub_args, genes_eligible_for_sub_args)    


def write_sub_arg_barcode_summary(
    df_arg, genes_eligible_for_sub_args, sub_arg_seqs_by_barcode, arg_names, 
    sub_arg_seqs_list, sub_arg_barcode_summary_tsv, extra_mle_info_sub_arg_tsv,
    p_match, p_none, p_error, alpha_prior, beta_prior,
    min_confidence, min_noise_reads, noise_cutoff_ratio):
    
    df_sub_arg = pd.DataFrame(df_arg)
    df_sub_arg[arg_names] = df_sub_arg[arg_names].astype(str)

    barcodes = df_arg.index
    mle_output_cols = ["Predicted_sub-ARG", "Confidence", "Contamination", "Total_#_of_ARG_reads", "Technical_noise_count"]
    cols_for_mle_df = []
    for gene in genes_eligible_for_sub_args:
        for mle_col_name in mle_output_cols:
            cols_for_mle_df.append(f"{gene}: {mle_col_name}")
    df_extra_mle_info = pd.DataFrame(index = barcodes, columns = cols_for_mle_df)
    df_extra_mle_info[:] = "-"

    for bc in barcodes:
        read_counts_all_args = sub_arg_seqs_by_barcode[bc]
        
        for arg in genes_eligible_for_sub_args:

            gene_num = arg_names.index(arg)

            if df_arg.loc[bc, arg] == 0:
                continue

            read_counts = {}
            for seq, count in read_counts_all_args.items():
                if int(seq.split("_")[1]) == gene_num:
                    read_counts[seq.split("_")[0]] = count

            final_sub_arg, total_reads, technical_noise_count, last_valid_confidence, bayesian_contamination_mean \
            = sub_arg_parse_and_analyze_perfect_corrected(
            read_counts, p_match, p_none, p_error, alpha_prior, beta_prior,
            min_confidence, min_noise_reads, noise_cutoff_ratio)

            if final_sub_arg is not None:
                df_sub_arg.loc[bc, arg] = final_sub_arg
                df_extra_mle_info.loc[bc, f"{arg}: Predicted_sub-ARG"] = final_sub_arg
                df_extra_mle_info.loc[bc, f"{arg}: Confidence"] = last_valid_confidence
                df_extra_mle_info.loc[bc, f"{arg}: Contamination"] = bayesian_contamination_mean
                df_extra_mle_info.loc[bc, f"{arg}: Total_#_of_ARG_reads"] = total_reads
                df_extra_mle_info.loc[bc, f"{arg}: Technical_noise_count"] = technical_noise_count
            else:
                df_sub_arg.loc[bc, arg] = f"{arg}_parent"

    with open(sub_arg_seqs_list, 'w') as f:

        f.write("Sub-ARG_Arbitrary_Name\tCell_count\tCore_sequence\n")

        # arbitrarily name each of the final sub-ARG sequences
        # replace sequences with their new names in the barcode summary
        # also record these names-to-sequences in a text file
        for arg in genes_eligible_for_sub_args:
            final_sub_arg_seqs = df_sub_arg.loc[ \
                (df_sub_arg[arg] != f"{arg}_parent") & \
                (df_sub_arg[arg] != "0") & \
                (df_sub_arg[arg] != "0.0"), \
                arg].to_list()

            ranked_final_sub_arg_seqs = Counter(final_sub_arg_seqs).most_common()

            for i, (seq, cell_count) in enumerate(ranked_final_sub_arg_seqs, start=1):
                sub_arg_name = f"{arg}_<{get_alpha_name(i)}>"
                f.write(f"{sub_arg_name}\t{cell_count}\t{seq}\n")
                df_sub_arg.loc[df_sub_arg[arg] == seq, arg] = sub_arg_name

                df_extra_mle_info.loc[df_extra_mle_info[f"{arg}: Predicted_sub-ARG"] == seq, \
                    f"{arg}: Predicted_sub-ARG"] = sub_arg_name

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
                                
    df_sub_arg.to_csv(sub_arg_barcode_summary_tsv, sep = "\t", index_label = "Barcode")
    df_extra_mle_info.to_csv(extra_mle_info_sub_arg_tsv, sep = "\t", index_label = "Barcode")

    return df_sub_arg

def modify_sub_arg_barcode_summary(
    df_sub_arg, final_sub_arg_names, filtered_sub_arg_barcode_summary_tsv,
    baseline_gene, is_baseline_gene_eligible_for_sub_args, genes_eligible_for_sub_args):

    # if baseline gene is not actually eligible for sub-ARGs,
    # modify its column to include only the parent ARG name, instead of sub-ARG names
    if not is_baseline_gene_eligible_for_sub_args:
        df_sub_arg.loc[df_sub_arg[baseline_gene] != "0", baseline_gene] = baseline_gene

    # for cells containing sub-ARGs that were filtered out, replace the sub-ARG name with the parent ARG name
    for arg in genes_eligible_for_sub_args:
        df_sub_arg.loc[~df_sub_arg[arg].isin(["0", f"{arg}_parent"] + final_sub_arg_names), arg] = f"{arg}_parent"

    df_sub_arg.to_csv(filtered_sub_arg_barcode_summary_tsv, sep = "\t", index_label = "Barcode")


def get_alpha_name(sub_arg_num: int) -> str:
    "Determine alpha-based name for sub-ARG, based on the number name equivalent"

    # note that input sub-ARG numbers are 1-indexed, not 0-indexed
    alphabet = "abcdefghijklmnopqrstuvwxyz".upper()
    assert sub_arg_num < 26*27, \
        f"Error: there should be less than 26*27 = 702 sub-ARGs per ARG, for the arbitrary sub-ARG naming system to work"

    starting_char = ""
    if sub_arg_num > 26:
        starting_char = alphabet[(sub_arg_num - 1)//26 - 1]

    alpha_name = starting_char + alphabet[sub_arg_num - 1]
    return alpha_name

def match_ids_to_genes_bc(
    b_with_ids: str, arg_packets: str, 
    filtered_barcode_summary_tsv: str) -> Dict[str, int]:
    """
    Goal: obtain information for all ARG read IDs related to the current filtered barcode summary.
    Create a dictionary with elements of the form "ID: gene|barcode"

    This dictionary will contain all IDs that are both 1) ARG reads, 
     and 2) ID corresponds to one of the barcodes in the input filtered_barcode_summary_tsv 
     (i.e. some barcode filtering has already been done - this ensures we avoid ARG read with sequences that are too short,
      as that would affect later processing, specifically the extract_core function).

    This dictionary will be used in downstream processing.
    """
    
    filtered_bcs = pd.read_csv(filtered_barcode_summary_tsv, sep="\t", index_col = "Barcode").index

    filtered_arg_ids_with_bs = []
    with open(b_with_ids, 'r') as b_with_ids_file:
        for line in b_with_ids_file:
            bc, ids = line.strip().split(": ")
            if bc in filtered_bcs:
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

    filtered_arg_ids = [id_with_b.split("|")[0] for id_with_b in filtered_arg_ids_with_bs]

    ids_to_genes_bc = {}
    with open(arg_packets, 'r') as arg_packet_file:
        for line in arg_packet_file:
            packet = json.loads(line.strip())
            id = packet["ID"]
            if id in filtered_arg_ids:
                bc = filtered_arg_ids_with_bs[filtered_arg_ids.index(id)].split("|")[1]
                gene_num = packet["gene"].index(1)
                ids_to_genes_bc[id] = f"{gene_num}|{bc}"
    return ids_to_genes_bc

def get_all_sub_arg_seqs(
    ids_to_genes: Dict, arg_names: List[str], fwd_fastq: str, rev_fastq: str, 
    nums_genes_eligible_for_sub_args: List[int]) -> List[Dict]:
    """
    Obtain all sub-ARG sequences present, for all ARGs that are eligible to have sub-ARGs.

    Returned value is a list of dictionaries, where each dictionary corresponds to one ARG.
    Each dictionary is of the form: (e.g. of TEM as being the parent ARG for this dictionary)
    {"Sub-ARG seq #1": ["ID#5|Barcode#30", "ID#8|Barcode#7"],   # corresponds to e.g. TEM_<1>
     "Sub-ARG seq #2": ["ID#10|Barcode#1"], ...}                # corresponds to e.g. TEM_<2>
    Note that all values chosen for this example are random - the goal above was just to show
    the sub_arg_seqs formatting.
    Each ID corresponding to a specific sub-ARG sequence are added to that sub-ARG's list
    (e.g. ID#5 and ID#8 would correspond to TEM_<1> in the example above.
    The ID names are from the fastq file. Barcodes are included with the IDs in the format of
    "ID|Barcode", to faciliate later downstream processing.
    Note also that "Sub-ARG seq #1" would be in the format outputted by the extract_core() function;
    that is, "ATCG...|TACC..." (with the specific sequences provided here being just an example).

    Return a list of dictionaries mapping each sub-ARG sequence to all its relevant IDs,
    while tracking parent ARGs (for the sub-ARGs) and barcodes (for the IDs).
    """

    sub_arg_seqs = [{} for _ in range(len(arg_names))] 

    with open_maybe_gzip(fwd_fastq, 'r') as fwd, open_maybe_gzip(rev_fastq, 'r') as rev:
        # read in the read_ID line
        f_line = fwd.readline()
        r_line = rev.readline()

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

            gene_num_and_bc = ids_to_genes.get(id_f)
            if gene_num_and_bc != None:
                gene_num, bc = gene_num_and_bc.split("|")
                gene_num = int(gene_num)

                # skip over ARGs that are not eligible to have sub-ARGs
                if gene_num not in nums_genes_eligible_for_sub_args:
                    i += 1 # ensure counter still increments, before "continue" keyword
                    continue

                seq_pair = extract_core(f_seq, r_seq)
                if seq_pair in sub_arg_seqs[gene_num].keys():
                    sub_arg_seqs[gene_num][seq_pair].append(f"{id_f}|{bc}")
                else:
                    sub_arg_seqs[gene_num][seq_pair] = [f"{id_f}|{bc}"]
            i += 1

    return sub_arg_seqs

def cluster_sub_arg_seqs(
    arg_names: List[str], sub_arg_seqs: List[Dict], nums_genes_eligible_for_sub_args: List[int], 
    max_mm_sub_arg: int, max_shift_sub_arg: int) -> List[Dict]:
    """
    Cluster the sub-ARG sequences, to create a whitelist of sub-ARG sequences for each ARG.
    Sub-ARG sequences are clustered based on the given max mismatch (should be zero) and max shift.
    Returns the modified sub_arg_seqs list (now called clustered_sub_arg_seqs)
    """

    clustered_sub_arg_seqs = [{} for _ in range(len(arg_names))]
    for gene_num in range(len(sub_arg_seqs)):

        # skip over ARGs that are not eligible to have sub-ARGs
        if gene_num not in nums_genes_eligible_for_sub_args:
            continue

        counts = {}
        for seq, ids in sub_arg_seqs[gene_num].items():
            counts[seq] = len(ids)
        sorted_sub_seqs = [k for k, v in Counter(counts).most_common()]

        for seq_s in sorted_sub_seqs:
            in_dominant_seqs = False
            for seq_d in clustered_sub_arg_seqs[gene_num]: # iterate through all dominant sub-ARG sequences identified so far for that ARG
                d = semi_global_distance(seq_s, seq_d, max_shift=max_shift_sub_arg)
                if d <= max_mm_sub_arg and within_cluster_boundary(clustered_sub_arg_seqs[gene_num], seq_s, max_mm_sub_arg, max_shift_sub_arg):
                    in_dominant_seqs = True
                    clustered_sub_arg_seqs[gene_num][seq_d].extend(sub_arg_seqs[gene_num][seq_s]) 
                    seq_added = True
                    break # similar sequence already in dominant sequences list
            if not in_dominant_seqs:
                clustered_sub_arg_seqs[gene_num][seq_s] = sub_arg_seqs[gene_num][seq_s]
                seq_added = True
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
    
    sub_arg_seqs_1_dict = {}
    total_num_sub_arg_seqs = 0

    for i in range(len(clustered_sub_arg_seqs)):
        
        # note that i is the gene number
        new_one_arg_dict = {}

        for seq, ids in clustered_sub_arg_seqs[i].items():
            new_one_arg_dict[f"{seq}_{i}"] = ids

        sub_arg_seqs_1_dict.update(new_one_arg_dict)
        total_num_sub_arg_seqs += len(new_one_arg_dict)        

    # check that no sub-ARG sequences were lost when updating dict
    assert len(sub_arg_seqs_1_dict) == total_num_sub_arg_seqs 

    return sub_arg_seqs_1_dict

def sort_sub_args_by_barcode(
    df_arg, sub_arg_seqs_1_dict: Dict) -> Dict:

    filtered_bcs = df_arg.index
    sub_arg_seqs_by_barcode = dict.fromkeys(filtered_bcs)
    for bc in sub_arg_seqs_by_barcode:
        sub_arg_seqs_by_barcode[bc] = {}

    for seq, ids_with_bs in sub_arg_seqs_1_dict.items():
        for id_and_b in ids_with_bs:
            bc = id_and_b.split("|")[1]
            if seq in sub_arg_seqs_by_barcode[bc]:
                sub_arg_seqs_by_barcode[bc][seq] += 1
            else:
                sub_arg_seqs_by_barcode[bc][seq] = 1

    return sub_arg_seqs_by_barcode

def get_genes_eligible_for_sub_args(
    primers_file: str) -> List[str]:
    """
    Extract information from primers file to form a list of all genes eligible to have sub-ARGs.
    The ARGs in the primers file are classified as either "single" or "family".
    Those classified as "family" are eligible to have sub-ARGs.
    Returns a list of eligible genes.
    """
    genes_eligible_for_sub_args = []
    with open(primers_file, 'r') as f:
        i = 0
        line = f.readline()
        while line != "":
            if i != 0 and i != 1: # skip over header line and 16s primer; only include the ARG primers
                if line.split(",")[3].strip().lower() == "family":
                    genes_eligible_for_sub_args.append(line.split(",")[0])
            line = f.readline()
            i += 1

    return genes_eligible_for_sub_args

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

def within_cluster_boundary(
    candidate_seqs: List[str], new_seq: str, 
    max_dist=0, max_shift_sub_arg=2) -> bool:
    """
    Check whether adding new_seq keeps the whole cluster boundary <= max_dist.
    """
    for old_seq in candidate_seqs:
        d = semi_global_distance(new_seq, old_seq, max_shift=max_shift_sub_arg)
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
    Infer one corrected taxonomy path and contamination estimate per barcode.

    The algorithm works in three broad stages:
    1. Use maximum likelihood to choose the dominant taxon at each rank.
    2. Stop going deeper when the best rank assignment is not confident enough.
    3. Treat small off-path taxa as technical noise and larger off-path taxa as
       real contamination before computing a Bayesian contamination estimate.
    """

    total_reads = sum(read_counts.values())

    # ---------------- 1. Stepwise adaptive MLE taxonomy path inference ----------------
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

    # If confidence is high enough, accept sub-ARG rank with the best candidate sub-ARG.
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
    parser.add_argument("--filtered_sub_arg_barcode_summary_tsv", type=str, default="reports/cell_amr_matrix.tsv")
    parser.add_argument("--baseline_gene", type=str, required=True)
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


    args = parser.parse_args()

    # 2026-08-10: Materialize temporary and final report directories before writing sub-ARG tables.
    # Reason: the cell-by-AMR matrix is a report while working tables remain under tmp.
    ensure_output_directories(
        args.sub_arg_barcode_summary_tsv, args.sub_arg_seqs_list,
        args.filtered_sub_arg_barcode_summary_tsv,
        args.filtered_stats_cells_per_sub_arg_tsv, args.extra_mle_info_sub_arg_tsv)
    
    # make sure input file paths exist
    if not os.path.exists(args.filtered_counts_summary_arg_tsv):
        print(f"❌ Error: input file not found: {args.filtered_counts_summary_arg_tsv}")
        return
    if not os.path.exists(args.r1_fastq):
        print(f"❌ Error: input file not found: {args.r1_fastq}")
        return
    if not os.path.exists(args.r2_fastq):
        print(f"❌ Error: input file not found: {args.r2_fastq}")
        return
    if not os.path.exists(args.primers_file):
        print(f"❌ Error: input file not found: {args.primers_file}")
        return
    if not os.path.exists(args.b_with_ids):
        print(f"❌ Error: input file not found: {args.b_with_ids}")
        return
    if not os.path.exists(args.arg_packets):
        print(f"❌ Error: input file not found: {args.arg_packets}")
        return


    create_sub_arg_barcode_summary(
        args.filtered_counts_summary_arg_tsv, 
        args.b_with_ids, args.arg_packets, 
        args.r1_fastq, args.r2_fastq,
        args.sub_arg_barcode_summary_tsv, args.primers_file, 
        args.sub_arg_seqs_list, 
        args.filtered_sub_arg_barcode_summary_tsv, 
        args.baseline_gene, args.filtered_stats_cells_per_sub_arg_tsv, 
        args.extra_mle_info_sub_arg_tsv,
        args.p_match, args.p_none, args.p_error, args.alpha_prior, args.beta_prior,
        args.min_confidence, args.min_noise_reads, args.noise_cutoff_ratio,
        args.alpha, args.max_shift_sub_arg, args.max_mm_sub_arg)

if __name__ == "__main__":
    main()
