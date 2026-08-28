#!/usr/bin/env python3
import sys
import re
import pandas as pd
from collections import defaultdict, Counter
from typing import List, Dict
from helper_functions import get_arg_names, ensure_output_directories
import os
import argparse
# 2026-08-28: Exit non-zero when a required stage input is missing.
# Reason: the public CLI uses subprocess check=True to stop failed stages.

# =====================================================================
# CORE BIOPHYSICAL & BIOINFORMATICS HYPERPARAMETERS
# =====================================================================
# Intragenomic variation tolerance: 16S rRNA genes within the same bacterial 
# genome can harbor 1-2 bp natural polymorphism due to multiple operon copies.
MAX_DISTANCE = 3              

# Maximum permitido shift (bp) during semi-global sequence alignment to handle 
# indexing/sequencing length variations.
MAX_SHIFT = 3                 

# Minimum sequencing depth (read count) required to evaluate a single-cell droplet.
MIN_READS = 5                 

# Bio-physical threshold for contamination: If a distinct out-of-bounds 
# (>2bp) secondary ASV accounts for >= 20% of the valid read pool, it indicates 
# a true doublet/co-encapsulation event rather than ambient background DNA soup.
MIXED_RATIO_THRESHOLD = 0.10  

# 1-based extraction coordinates to accurately truncate and merge R1/R2 reads 
# while safely skipping the 20bp cell barcode sequence at the start of R2.
R1_START, R1_END = 30, 120
R2_START, R2_END = 70, 120 


def semi_global_distance(a, b, max_shift=2):
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

def within_cluster_boundary(candidate_seqs, new_seq, max_dist=2):
    """
    Check whether adding new_seq keeps the whole cluster boundary <= max_dist.
    """
    for old_seq in candidate_seqs:
        d = paired_semi_global_distance(new_seq, old_seq, max_shift=MAX_SHIFT)
        if d > max_dist:
            return False, d
    return True, max_dist


def read_final_barcodes(barcode_summary_tsv_filename: str) -> List[str]:
    """Returns all barcodes from the filtered barcode summary tsv table"""
    
    df = pd.read_csv(barcode_summary_tsv_filename, sep="\t", index_col = "Barcode")
    return df.index.to_list()


def read_b_with_ids(b_with_ids_filename: str) -> Dict[str, List[str]]:
    """Maps each cell barcode to its corresponding 16s sequencing Read IDs."""

    bc_to_ids = {}

    with open(b_with_ids_filename, 'r') as b_with_ids_file:
        for line in b_with_ids_file:
            # parse barcode line
            bc, ids = line.strip("\n").split(": ")
            _16s_ids = ids.split(" | ")[0].split(", ")

            # 2026-08-27: Normalize and remove empty read IDs in one pass.
            # Reason: the previous exact-string removals missed tabs and other whitespace-only IDs.
            _16s_ids = [read_id.strip() for read_id in _16s_ids if read_id.strip()]

            # add to bc_to_ids dict
            bc_to_ids[bc] = _16s_ids
            
    return bc_to_ids


def read_16s_revised(fwd_16s_fastq: str, rev_16s_fastq: str) -> Dict[str, List[str]]:
    """Read in the sequences for all 16s reads"""
    _16s_reads = {}

    with open(fwd_16s_fastq, 'r') as fwd_file, open(rev_16s_fastq, 'r') as rev_file:
        # reading in and parsing fastq files
        f_line = fwd_file.readline()
        r_line = rev_file.readline()
        # 2026-08-27: Track the current FASTQ record for mismatch diagnostics.
        # Reason: the previous error message referenced undefined variable i.
        line_num = 1

        while (f_line != "") and (r_line != ""):
            id_f = f_line.strip().split(" ")[0].strip("@")
            id_r = r_line.strip().split(" ")[0].strip("@")
            assert id_f == id_r, f"R1/R2 read ID mismatch near FASTQ line {line_num}"
            f_seq = fwd_file.readline().strip()
            r_seq = rev_file.readline().strip()
            
            # main functionality: add to _16s_reads dict
            _16s_reads[id_f] = [f_seq, r_seq]

            # prepare for next read, accounting for 4-line fastq format
            for _ in range(3):
                f_line = fwd_file.readline()
                r_line = rev_file.readline()
            line_num += 4

    return _16s_reads
    

def extract_core(r1, r2):
    """Slices and concatenates the predefined hypervariable core regions from R1 and R2."""
    if len(r1) < R1_END or len(r2) < R2_END:
        return None
    return r1[R1_START:R1_END] + "|" + r2[R2_START:R2_END]


# =====================================================================
# CORE ALGORITHM: GRADIENT THREE-TIER CLASSIFICATION ENGINE
# =====================================================================
def summarize_barcode(core_counter):
    """
    Performs single-cell de-multiplexing and quality control using a centered 
    dominant ASV topology. Handles intragenomic heterogeneity (<=2bp) and filters 
    ambient contamination via a progressive bio-physical ratio matrix.
    """
    reads_used = sum(core_counter.values())

    if reads_used == 0:
        return {
            "reads_used": 0, "raw_unique": 0, "dominant_raw_count": 0,
            "coexisting_2bp_reads": 0, "final_count": 0, "max_observed_dist": 0,
            "secondary_reads": 0, "status": "no_usable_reads", "consensus_seq": "NA"
        }

    # 1. Lock onto the Absolute Dominant Sequence (The true biological center of the cell)
    dominant_seq, raw_dominant_count = core_counter.most_common(1)[0]
    
    # 2. Dynamic low-frequency noise gate (Filters 5% baseline sequencer cross-talk/errors)
    noise_threshold = max(2, int(reads_used * 0.05))

    final_count = raw_dominant_count
    coexisting_2bp_reads = 0
    max_observed_dist = 0
    cluster_seqs = [dominant_seq]
    
    # Tracks reads of distinct out-of-bounds (>2bp) secondary alleles that clear the noise gate
    total_unauthorized_secondary_reads = 0

    # 3. Topological Clustering and Sorting Matrix
    for seq, count in core_counter.most_common():
        if seq == dominant_seq:
            continue
            
        # Compute sliding distance to the single dominant hub
        d = paired_semi_global_distance(seq, dominant_seq, max_shift=MAX_SHIFT)
        
        if d <= MAX_DISTANCE:
            ok, boundary_d = within_cluster_boundary(
                cluster_seqs,
                seq,
                max_dist=MAX_DISTANCE
            )   

            if ok:
                final_count += count
                coexisting_2bp_reads += count
                cluster_seqs.append(seq)

                if d > max_observed_dist:
                    max_observed_dist = d
            else:
                if count >= noise_threshold:
                    total_unauthorized_secondary_reads += count
                    if boundary_d > max_observed_dist:
                        max_observed_dist = boundary_d
        else:
            if count >= noise_threshold:
                total_unauthorized_secondary_reads += count
                if d > max_observed_dist:
                    max_observed_dist = d

    # 4. Bio-Physical Ratio Classification (Three-Tier Resolution)
    total_valid_pool = final_count + total_unauthorized_secondary_reads
    secondary_ratio = total_unauthorized_secondary_reads / total_valid_pool if total_valid_pool > 0 else 0.0

    if reads_used < MIN_READS:
        status = "low_depth"
    elif total_unauthorized_secondary_reads == 0:
        # Perfectly pristine single cell with no external interference
        status = "single_ASV"
    elif secondary_ratio >= MIXED_RATIO_THRESHOLD:
        # High-abundance alien strain detected (>= 20%). Conclusive evidence of microfluidic doublet/multiplet.
        status = "mixed_ASV"
    else:
        # [THE MIDDLE GROUND]: Distinct ASV exists but represents < 20% of the droplet.
        # Categorized as an ambient DNA soup contamination or a rare structural copy number outlier.
        # Retained as a functional single cell while flagging the technical imperfection.
        status = "corrupted_single_ASV"

    final_fraction = final_count / reads_used if reads_used > 0 else 0.0

    return {
        "reads_used": reads_used,
        "raw_unique": len(core_counter),
        "dominant_raw_count": raw_dominant_count,
        "coexisting_2bp_reads": coexisting_2bp_reads,
        "final_count": final_count,
        "final_fraction": final_fraction,
        "max_observed_dist": max_observed_dist,
        "secondary_reads": total_unauthorized_secondary_reads,
        "status": status,
        "consensus_seq": dominant_seq
    }



# =====================================================================
# DOWNSTREAM REPORT DATA EXPANSION AND GENERATION
# =====================================================================
def conduct_asv_typing(
    barcode_summary_tsv_filename: str, b_with_ids_filename: str,
    fwd_16s_fastq: str, rev_16s_fastq: str,
    asv_barcode_summary_tsv_filename: str, global_asv_tsv_filename: str,
    primers_file: str, filter_corrupted: bool = False
    ):
    """Run the full ASV typing algorithm"""

    # read in input files
    final_barcodes = read_final_barcodes(barcode_summary_tsv_filename)
    bc_to_ids = read_b_with_ids(b_with_ids_filename)
    _16s_reads = read_16s_revised(fwd_16s_fastq, rev_16s_fastq)

    # run ASV typing for each cell
    barcode_summary = {}
    global_counter = Counter()
    for bc in final_barcodes:
        core_counter = Counter()
        for rid in bc_to_ids.get(bc, []):
            core = extract_core(_16s_reads[rid][0], _16s_reads[rid][1])
            if core is None:
                continue
            core_counter[core] += 1

        s = summarize_barcode(core_counter)
        barcode_summary[bc] = s

        # Both pristine and corrupted single cells contribute their dominant ASV to the global registry
        if s["status"] in ["single_ASV", "corrupted_single_ASV"]:
            global_counter[s["consensus_seq"]] += 1


    # 2026-08-28: Emit standardized public ASV cluster identifiers.
    # Reason: distinguish generated sequence clusters from biological variant names.
    # write global ASV output file, which contains all the ASV sequences matched to their names (e.g. ASV_1)
    seq_to_asv = {}
    with open(global_asv_tsv_filename, "w") as out:
        out.write("Core_ASV_ID\tcell_count\tcore_sequence\n")
        for i, (seq, cell_count) in enumerate(global_counter.most_common(), start=1):
            asv_id = f"ASV_{i}"
            seq_to_asv[seq] = asv_id
            out.write(f"{asv_id}\t{cell_count}\t{seq}\n")

    # write barcode summary, but with new ASV columns appended
    # also filter out cells classified as having mixed ASVs
    write_ASV_barcode_summary(
        barcode_summary_tsv_filename, final_barcodes, \
        barcode_summary, seq_to_asv, \
        asv_barcode_summary_tsv_filename, primers_file, filter_corrupted)


def filter_ASV(df_combined, filter_corrupted: bool = False):
    """Filter out cells classified as having mixed ASVs, and optionally low-confidence single ASVs"""
    
    original_num_barcodes = len(df_combined)

    extra_info = ""
    if filter_corrupted: 
        extra_info = " or corrupted single ASV"
    print("Filter out barcodes identified as having a mixed ASV" + extra_info +  "...")

    # 2026-08-10: Build one removal mask before changing the DataFrame.
    # Reason: deleting a row during iteration can make the next .loc lookup fail.
    statuses_to_remove = {"mixed_ASV", "low_depth"}
    if filter_corrupted:
        statuses_to_remove.add("corrupted_single_ASV")
    df_combined.drop(
        index=df_combined.index[df_combined["Status"].isin(statuses_to_remove)],
        inplace=True,
    )

    print("Pre-filtering # of barcodes:", original_num_barcodes)
    print("  # of barcodes filtered out:", original_num_barcodes-len(df_combined))
    print("  # of barcodes remaining:", len(df_combined))


def write_ASV_barcode_summary(filtered_barcode_summary_tsv_filename: str, \
    final_barcodes: List[str], barcode_summary: Dict, seq_to_asv: Dict, \
    asv_barcode_summary_tsv_filename: str, primers_file: str, filter_corrupted: bool = False):
    """
    Write revised barcode summary - specifically, add new columns for per-cell ASV information, 
    and filter cells based on their ASV status.
    """
    
    df_original = pd.read_csv(filtered_barcode_summary_tsv_filename, sep="\t", index_col = "Barcode")
    df_asv = pd.DataFrame(columns = ["Reads_used_for_ASV", "Raw unique_core_sequences", \
    "Dominant_raw_read_count", "Coexisting_2bp_reads", "Unauthorized_secondary_reads", \
    "Final_cell_asv_reads", "Max_internal_distance", "Assigned_core_asv", "Status"])

    arg_names = get_arg_names(primers_file)
    df_original[arg_names] = df_original[arg_names].astype("Int64")

    # write ASV-specific columns of barcode summary
    i = 0
    for bc in final_barcodes:
        s = barcode_summary[bc]
        asv_id = seq_to_asv.get(s["consensus_seq"], "NA")
        df_asv.loc[bc] = [s['reads_used'], s['raw_unique'], s['dominant_raw_count'], \
        s['coexisting_2bp_reads'], s['secondary_reads'], s['final_count'], \
        s['max_observed_dist'], asv_id, s['status']]

        i += 1
        if i % 1000 == 0:
            print(f"Processed {i//1000},000 barcodes")

    # re-arrange order of columns in barcode summary to be: MLE info | ASV info | ARG info
    mle_info_cols = ["Predicted taxonomy", "Confidence", "Contamination", "Total # of 16s reads", "Technical noise count"]
    df_original_mle_info = df_original[mle_info_cols]
    df_original_args = pd.DataFrame(df_original)
    for col in mle_info_cols:
        df_original_args.drop(col, axis = 1, inplace = True)

    df_combined = pd.concat([df_original_mle_info, df_asv, df_original_args], axis=1)
    filter_ASV(df_combined, filter_corrupted) # filter out cells with mixed and optionally low-confidence single ASVs

    df_combined.to_csv(asv_barcode_summary_tsv_filename, sep = "\t", index_label = "Barcode")


# 2026-08-10: Parse command-line booleans from explicit true/false strings.
# Reason: bool("False") is True in Python and silently enables the wrong ASV filter.
def parse_bool(value: str) -> bool:
    normalized = value.strip().lower()
    if normalized in {"true", "1", "yes", "y"}:
        return True
    if normalized in {"false", "0", "no", "n"}:
        return False
    raise argparse.ArgumentTypeError("expected true/false")


def main():
    parser = argparse.ArgumentParser()

    # take input parameters
    parser.add_argument("--barcode_summary_tsv", type=str, required=True)
    parser.add_argument("--b_with_ids", type=str, required=True)
    # 2026-08-10: Expose ASV paired reads as R1/R2 in the public CLI.
    # Reason: sequencing inputs are named R1 and R2, not forward/reverse.
    parser.add_argument("--r1_16s_fastq", type=str, required=True)
    parser.add_argument("--r2_16s_fastq", type=str, required=True)
    # 2026-08-10: Route sequence-variant tables to tmp by default.
    # Reason: these tables support taxonomic assignment but are not the primary result.
    parser.add_argument("--asv_barcode_summary_tsv", type=str, default = "tmp/asv_barcode_summary.tsv")
    parser.add_argument("--global_asv_tsv", type=str, default = "tmp/global_asv.tsv")
    # 2026-08-10: Use parse_bool for predictable CLI behavior.
    # Reason: argparse type=bool treats every non-empty string as True.
    parser.add_argument("--filter_corrupted", type=parse_bool, default=False)
    parser.add_argument("--primers_file", type=str, required=True)

    args = parser.parse_args()

    # 2026-08-10: Materialize the temporary output directory before writing ASV tables.
    # Reason: the default tmp paths must work in a new result directory.
    ensure_output_directories(args.asv_barcode_summary_tsv, args.global_asv_tsv)
    
    # make sure input file paths exist
    if not os.path.exists(args.barcode_summary_tsv):
        print(f"❌ Error: input file not found: {args.barcode_summary_tsv}")
        sys.exit(1)
    if not os.path.exists(args.b_with_ids):
        print(f"❌ Error: input file not found: {args.b_with_ids}")
        sys.exit(1)
    if not os.path.exists(args.r1_16s_fastq):
        print(f"❌ Error: input file not found: {args.r1_16s_fastq}")
        sys.exit(1)
    if not os.path.exists(args.r2_16s_fastq):
        print(f"❌ Error: input file not found: {args.r2_16s_fastq}")
        sys.exit(1)
    if not os.path.exists(args.primers_file):
        print(f"❌ Error: input file not found: {args.primers_file}")
        sys.exit(1)

    conduct_asv_typing(
        args.barcode_summary_tsv, args.b_with_ids,
        args.r1_16s_fastq, args.r2_16s_fastq,
        args.asv_barcode_summary_tsv, args.global_asv_tsv, 
        args.primers_file, args.filter_corrupted
    )

if __name__ == "__main__":
    main()
