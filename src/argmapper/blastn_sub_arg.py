import subprocess
import os
import argparse
import pandas as pd
MAX_HITS = 500

def make_summary_table(
    sub_arg_seqs_list, blastn_sub_arg_tsv, query_fasta, input_fasta, db,
    final_barcode_summary_tsv, first_gene_column_num):
        
    make_blast_db(input_fasta, db) # make the BLAST database

    final_gene_names = get_final_gene_names(final_barcode_summary_tsv, first_gene_column_num)

    header = ["sub-ARG_name", "subject",
              "read1_start", "read1_end", "read1_percent_identity",
              "read2_start", "read2_end", "read2_percent_identity",
              "read1_seq", "read2_seq"]

    with open(sub_arg_seqs_list, 'r') as input, \
    open(blastn_sub_arg_tsv, 'w') as output:
        output.write("\t".join(header) + "\n")
        
        i = 0
        for line in input:
            i += 1
            if line.split(": ")[0] in final_gene_names:
                name, seqs = line.strip().split(": ")
                f_seq, r_seq = seqs.split("|")
                rows = get_blastn_rows(name, f_seq, r_seq, query_fasta, db)
                for row in rows:
                    output.write("\t".join(row) + "\n")

def get_final_gene_names(final_barcode_summary_tsv, first_gene_column_num):

    df_summary = pd.read_csv(final_barcode_summary_tsv, sep="\t", index_col = "Barcode")
    all_col_names = df_summary.columns.to_list()
    final_gene_names = all_col_names[first_gene_column_num:] 
    return final_gene_names

def make_blast_db(input_fasta, db):
    make_blast_db_command = [
        "makeblastdb",
        "-in", str(input_fasta),
        "-dbtype", "nucl",
        "-out", str(db)
    ]
    subprocess.check_output(make_blast_db_command, text=True)

def get_blastn_rows(sub_arg_name, read1_seq, read2_seq, query_fasta, db):
    
    # 1. Write the two reads to a FASTA file for BLAST
    with open(query_fasta, 'w') as f:
        f.write(f">read1\n{read1_seq}\n>read2\n{read2_seq}\n")

    # 2. Run BLAST 
    blast_columns = [
        "qseqid", "sseqid", "pident", "qstart", "qend", "qlen",
        "sstart", "send", "gapopen", "stitle",
    ]

    blast_command = [
        "blastn",
        "-query", str(query_fasta),
        "-db", db,
        "-max_target_seqs", str(MAX_HITS),
        "-outfmt", "6 " + " ".join(blast_columns),
    ]

    print(sub_arg_name)

    blast_output = subprocess.check_output(blast_command, text=True)
    print("blast output:\n\n", blast_output)

    with open(query_fasta, 'r') as f:
        print("query fasta:", f.read())

    # 3. Read BLAST hits and keep the best hit
    hits_by_subject = {}
    query_order = []

    for line in blast_output.splitlines():
        values = line.split("	")
        hit = dict(zip(blast_columns, values))

        query_name = hit["qseqid"]
        subject_name = hit["sseqid"]

        if query_name not in query_order:
            query_order.append(query_name)

        hits_by_subject.setdefault(subject_name, {})
        hits_by_subject[subject_name].setdefault(query_name, hit)

    # 4. Full length and 100%.
    rows = []
    read1_name, read2_name = query_order[:2]

    for subject_hits in hits_by_subject.values():
        if read1_name not in subject_hits or read2_name not in subject_hits:
            continue

        read1 = subject_hits[read1_name]
        read2 = subject_hits[read2_name]

        read1_aligned_length = abs(int(read1["qend"]) - int(read1["qstart"])) + 1
        read2_aligned_length = abs(int(read2["qend"]) - int(read2["qstart"])) + 1

        read1_is_exact = (
            read1_aligned_length == int(read1["qlen"])
            and float(read1["pident"]) == 100.0
            and int(read1["gapopen"]) == 0
        )
        read2_is_exact = (
            read2_aligned_length == int(read2["qlen"])
            and float(read2["pident"]) == 100.0
            and int(read2["gapopen"]) == 0
        )

        # only keep hits with exact matches on both fwd & rev reads
        if not (read1_is_exact and read2_is_exact):
            continue

        read1_start = str(min(int(read1["sstart"]), int(read1["send"])))
        read1_end = str(max(int(read1["sstart"]), int(read1["send"])))
        read2_start = str(min(int(read2["sstart"]), int(read2["send"])))
        read2_end = str(max(int(read2["sstart"]), int(read2["send"])))

        rows.append([
            sub_arg_name, read1["stitle"],
            read1_start, read1_end, "100%",
            read2_start, read2_end, "100%",
            read1_seq, read2_seq
        ])

    # for each sub-ARG, sort hit matches alphabetically by species name
    rows.sort(key=lambda row: row[1]) 

    return rows


def main():
    parser = argparse.ArgumentParser()

    # take input parameters
    parser.add_argument("--sub_arg_seqs_list", type=str, required=True)
    parser.add_argument("--blastn_sub_arg_tsv", type=str, default="blastn_sub_arg.tsv")
    parser.add_argument("--query_fasta", type=str, default="query_reads.fa")
    parser.add_argument("--input_fasta", type=str, required=True)
    parser.add_argument("--db", type=str, default="blastn_db")
    parser.add_argument("--final_barcode_summary_tsv", type=str, required=True)
    parser.add_argument("--first_gene_column_num", type=int, required=True)
    
    args = parser.parse_args()
    
    # make sure input file paths exist
    if not os.path.exists(args.sub_arg_seqs_list):
        print(f"❌ Error: input file not found: {args.sub_arg_seqs_list}")
        return
    if not os.path.exists(args.input_fasta):
        print(f"❌ Error: input file not found: {args.input_fasta}")
        return
    if not os.path.exists(args.final_barcode_summary_tsv):
        print(f"❌ Error: input file not found: {args.final_barcode_summary_tsv}")
        return

    make_summary_table(
    args.sub_arg_seqs_list, args.blastn_sub_arg_tsv, args.query_fasta, args.input_fasta, args.db,
    args.final_barcode_summary_tsv, args.first_gene_column_num)

if __name__ == "__main__":
    main()
