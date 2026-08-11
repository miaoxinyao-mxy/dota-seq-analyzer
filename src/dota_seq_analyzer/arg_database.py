import pandas as pd
import json
from collections import defaultdict, Counter
from typing import List, Dict, Tuple
from filter_sub_args import run_sub_arg_denoising_pipeline
from helper_functions import open_maybe_gzip, get_arg_names, ensure_output_directories
import os
import argparse

# 1-based extraction coordinates to accurately truncate and merge R1/R2 reads 
# while safely skipping the 20bp cell barcode sequence at the start of R2, as well as the ~20bp overlap sequence.
R1_START, R1_END = 30, 110
R2_START, R2_END = 70, 110 


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

def count_cells_per_sub_arg(
    sub_arg_barcode_summary_tsv: str, sub_arg_names: List[str], 
    stats_cells_per_sub_arg_tsv: str):

    df = pd.read_csv(sub_arg_barcode_summary_tsv, sep="\t", index_col = "Barcode")
    df_stats = pd.DataFrame(index = sub_arg_names, columns = ["Original_#_of_Cells"])
    df_stats[:] = 0

    for idx in df.index.to_list():
        for sub_arg in sub_arg_names:            
            k = df.loc[idx, sub_arg] # raw number of reads of the specific arg in this cell
            if k > 0:
                # record that this gene was originally positive for this cell
                df_stats.loc[sub_arg, "Original_#_of_Cells"] += 1 

    df_stats.to_csv(stats_cells_per_sub_arg_tsv, sep = "\t", index_label = "Sub-ARG")

def merge_into_one_dict(
    clustered_sub_arg_seqs: List[Dict]) -> Dict:
    
    sub_arg_seqs_1_dict = {}
    total_num_sub_arg_seqs = 0

    for arg_dict in clustered_sub_arg_seqs:
        sub_arg_seqs_1_dict.update(arg_dict)
        total_num_sub_arg_seqs += len(arg_dict)

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

def write_sub_arg_seqs_list_file(
    sub_arg_seqs_list: str, arg_names: List[str], 
    clustered_sub_arg_seqs: List[Dict]) -> List[str]:

    all_sub_args = []
    with open(sub_arg_seqs_list, 'w') as f:
        for i in range(len(arg_names)):
            arg_name = arg_names[i]
            f.write(arg_name + "\n\n")
            sub_arg_seqs = clustered_sub_arg_seqs[i]
            j = 1
            for seq in sub_arg_seqs:
                sub_arg_name = f"{arg_name}_<{j}>"
                all_sub_args.append(f"{sub_arg_name}:{seq}")
                f.write(f"{sub_arg_name}: {seq}\n")
                j += 1
            f.write("---------------------------------------------------------------------------------\n\n")
    return all_sub_args

def write_sub_arg_barcode_summary(
    df_arg, arg_names: List[str], all_sub_args: List[str], sub_arg_seqs_by_barcode: Dict, 
    sub_arg_barcode_summary_tsv: str, genes_eligible_for_sub_args: List[str]) -> Tuple:
    
    df_sub_arg = pd.DataFrame(df_arg)

    # drop names of ARGs eligible for sub-ARGs; these will be replaced with the sub-ARG names
    for arg_name in arg_names:
        if arg_name in genes_eligible_for_sub_args:
            df_sub_arg.drop(columns = arg_name, inplace = True)
    
    sub_arg_names = [x.split(":")[0] for x in all_sub_args] 
    df_new_cols = pd.DataFrame(index = df_sub_arg.index, columns = sub_arg_names)
    df_new_cols[:] = 0

    for bc in df_new_cols.index:
        seqs_for_b = sub_arg_seqs_by_barcode[bc]
        sub_arg_counts = [0 for _ in range(len(all_sub_args))]
        for i in range(len(all_sub_args)):
            if df_arg.loc[bc, all_sub_args[i].split("_")[0]] == 0:
                # if ARG count is zero or already filtered out, 
                # then sub-ARG count for that ARG would also be zero or filtered out
                continue
            sub_arg_seq = all_sub_args[i].split(":")[1]
            if sub_arg_seq in seqs_for_b:
                sub_arg_counts[i] += seqs_for_b[sub_arg_seq]
        df_new_cols.loc[bc] = sub_arg_counts
        

    df_sub_arg = pd.concat([df_sub_arg, df_new_cols], axis=1)
    df_sub_arg.to_csv(sub_arg_barcode_summary_tsv, sep = "\t", index_label = "Barcode")

    return df_sub_arg, sub_arg_names

def get_final_gene_names(
    arg_names: List[str], baseline_gene: str, genes_eligible_for_sub_args: List[str], 
    final_sub_arg_names: List[str], is_baseline_gene_eligible_for_sub_args: bool) -> List[str]:
    """Compile ordered list of final ARG & sub-ARG names"""

    final_gene_names = [] # names of all final sub-ARGs and ARGs
                          # ordered the same as original primers file (and thus arg_names list)
                          # (note that parent ARG names are only used for ARGs that are not eligible for sub-ARGs)
    for arg in arg_names:
        if arg in genes_eligible_for_sub_args:
            if arg == baseline_gene and not is_baseline_gene_eligible_for_sub_args:
                # baseline gene was included in the genes eligible for sub-ARGs
                # however, this was only so that we'd have a baseline gene to use in sub-ARG filtering
                # since this baseline gene is not actually eligible for sub-ARGs though,
                # make sure it's treated as such - i.e. treat as parent ARG, no sub-ARGs
                final_gene_names.append(arg)
                continue
            current_gene_sub_arg_names = [x for x in final_sub_arg_names if x.split("_")[0] == arg]
            current_gene_sub_arg_names.append(f"{arg}_<untyped>")
            final_gene_names.extend(current_gene_sub_arg_names)
        else:
            final_gene_names.append(arg)

    return final_gene_names

def write_filtered_sub_arg_barcode_summary(
    df_sub_arg, genes_eligible_for_sub_args: List[str], arg_names: List[str], sub_arg_names: List[str], 
    final_sub_arg_names: List[str], filtered_sub_arg_barcode_summary_tsv: str, baseline_gene: str, 
    is_baseline_gene_eligible_for_sub_args: bool, df_arg):

    final_gene_names = get_final_gene_names(arg_names, baseline_gene, genes_eligible_for_sub_args, 
                       final_sub_arg_names, is_baseline_gene_eligible_for_sub_args)    

    # Compile new final columns for ARGs & sub-ARGs
    df_new_cols = pd.DataFrame(index = df_sub_arg.index, columns = final_gene_names)

    for gene in final_gene_names:
        if gene == baseline_gene and not is_baseline_gene_eligible_for_sub_args:
            # need to use df_arg instead of df_sub_arg, because baseline gene was previously treated as a sub-ARG in df_sub_arg
            df_new_cols[gene] = df_arg[gene] 
            pass
        elif gene in df_sub_arg.columns:
            df_new_cols[gene] = df_sub_arg[gene]
        else:
            # account for sub-ARGs in the form of "arg_<untyped>"
            arg = gene.split("_")[0]
            removed_sub_args_for_gene = [x for x in df_sub_arg.columns if x.split("_")[0] == arg and x not in final_gene_names]
            df_temp = df_sub_arg[removed_sub_args_for_gene]
            df_new_cols[gene] = df_temp.sum(axis=1) 

    # Drop ARG and sub-ARG columns from original sub-ARG barcode summary
    original_gene_names = sub_arg_names + arg_names
    for gene_name in original_gene_names:
        if gene_name in df_sub_arg.columns:
            df_sub_arg.drop(columns = gene_name, inplace = True) 

    # Create the final barcode summary
    df_sub_arg = pd.concat([df_sub_arg, df_new_cols], axis=1)
    df_sub_arg.to_csv(filtered_sub_arg_barcode_summary_tsv, sep = "\t", index_label = "Barcode")


def create_sub_arg_barcode_summary(
    filtered_counts_summary_arg_tsv: str, 
    b_with_ids: str, arg_packets: str, 
    fwd_fastq: str, rev_fastq: str,
    sub_arg_barcode_summary_tsv: str, primers_file: str, 
    sub_arg_seqs_list: str, stats_cells_per_sub_arg_tsv: str, 
    filtered_sub_arg_barcode_summary_tsv: str, 
    baseline_gene: str, filtered_stats_cells_per_sub_arg_tsv: str, 
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
            
    # Write new barcode summary file, with counts for sub-ARGs instead of ARGs

    # 1) Determine the sub-ARG column names, and write file containing sequences for each sub-ARG name
    print("C. 1) Name each sub-ARG, and write the sub-ARG names to sequences in a file...")
    all_sub_args = write_sub_arg_seqs_list_file(sub_arg_seqs_list, arg_names, clustered_sub_arg_seqs)
    
    # 2) Write new dataframe barcode summary with sub-ARGs
    print("C. 2) Write barcode summary with columns as sub-ARGs (where applicable) instead of ARGs...")
    df_sub_arg, sub_arg_names = write_sub_arg_barcode_summary(
        df_arg, arg_names, all_sub_args, sub_arg_seqs_by_barcode, sub_arg_barcode_summary_tsv, genes_eligible_for_sub_args)

    # 3) Filter sub-ARGs in the barcode summary

    #  a) For each sub-ARG, determine # of cells containing at least 1 read of that sub-ARG
    print("C. 3)a) Determine # of cells containing at least 1 read of the given sub-ARG, for each sub-ARG...")
    count_cells_per_sub_arg(sub_arg_barcode_summary_tsv, sub_arg_names, stats_cells_per_sub_arg_tsv)

    #  b) Filter sub-ARGs, and obtain final names list of sub-ARGs to be used 
    print("C. 3)b) Filter sub-ARGs, and obtain final names list of sub-ARGs to be used...")
    final_sub_arg_names = run_sub_arg_denoising_pipeline(stats_cells_per_sub_arg_tsv, 
                        baseline_gene, alpha, filtered_stats_cells_per_sub_arg_tsv)

    #  c) Compile list of final ARG & sub-ARG names, and modify the barcode summary appropriately
    print("C. 3)c) Modify barcode summary based on filtered sub-ARGs...")
    print()
    write_filtered_sub_arg_barcode_summary(
        df_sub_arg, genes_eligible_for_sub_args, arg_names, sub_arg_names, 
        final_sub_arg_names, filtered_sub_arg_barcode_summary_tsv, baseline_gene, is_baseline_gene_eligible_for_sub_args, df_arg)


def main():
    parser = argparse.ArgumentParser()
  
    # take input parameters
    # 2026-08-10: Route sub-ARG working tables to tmp by default.
    # Reason: only the exported JSONL and selected reports are final products.
    parser.add_argument("--filtered_counts_summary_arg_tsv", type=str, default="tmp/filtered_counts_summary_arg.tsv")
    parser.add_argument("--b_with_ids", type=str, required=True)
    parser.add_argument("--arg_packets", type=str, required=True)
    # 2026-08-10: Expose ARG database inputs as R1/R2 in the public CLI.
    # Reason: sequencing inputs are named R1 and R2, not forward/reverse.
    parser.add_argument("--r1_fastq", type=str, required=True)
    parser.add_argument("--r2_fastq", type=str, required=True)
    parser.add_argument("--sub_arg_barcode_summary_tsv", type=str, default="tmp/sub_arg_barcode_summary.tsv")
    parser.add_argument("--primers_file", type=str, required=True)
    parser.add_argument("--sub_arg_seqs_list", type=str, default="tmp/sub_arg_seqs_list.txt")
    parser.add_argument("--stats_cells_per_sub_arg_tsv", type=str, default="tmp/stats_cells_per_sub_arg.tsv")
    # 2026-08-10: Use the generic cell-by-target report name.
    # Reason: the public workflow is no longer limited to AMR panels.
    parser.add_argument("--filtered_sub_arg_barcode_summary_tsv", type=str, default="reports/cell_target_matrix.tsv")
    parser.add_argument("--baseline_gene", type=str, required=True)
    parser.add_argument("--filtered_stats_cells_per_sub_arg_tsv", type=str, default="tmp/filtered_stats_cells_per_sub_arg.tsv")
    parser.add_argument("--alpha", type=float, default=0.05)
    parser.add_argument("--max_shift_sub_arg", type=int, default=2)
    parser.add_argument("--max_mm_sub_arg", type=int, default=0)

    args = parser.parse_args()

    # 2026-08-10: Materialize temporary and final report directories before writing sub-ARG tables.
    # Reason: the cell-by-AMR matrix is a report while working tables remain under tmp.
    ensure_output_directories(
        args.sub_arg_barcode_summary_tsv, args.sub_arg_seqs_list,
        args.stats_cells_per_sub_arg_tsv,
        args.filtered_sub_arg_barcode_summary_tsv,
        args.filtered_stats_cells_per_sub_arg_tsv)
    
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
        args.sub_arg_seqs_list, args.stats_cells_per_sub_arg_tsv, 
        args.filtered_sub_arg_barcode_summary_tsv, 
        args.baseline_gene, args.filtered_stats_cells_per_sub_arg_tsv, 
        args.alpha, args.max_shift_sub_arg, args.max_mm_sub_arg)


if __name__ == "__main__":
    main()
