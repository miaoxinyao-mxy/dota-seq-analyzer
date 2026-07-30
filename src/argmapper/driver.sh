#!/bin/bash
# Exit immediately if any command exits with a non-zero status
set -e

# ---------------------------------------------
# Define parameter values and filenames
# ---------------------------------------------

# Input paths
FWD_FASTQ="/home/julianna/workspace/pipeline_input_files/fastq_data/forward.fastq"
REV_FASTQ="/home/julianna/workspace/pipeline_input_files/fastq_data/reverse.fastq"
PRIMER_CSV="/home/julianna/workspace/pipeline_input_files/primer_dotaseq.csv"
DB_SUB_ARG_FASTA="/home/julianna/workspace/pipeline_input_files/db_sub_arg.fasta"

# Input yes/no to optional pipeline steps
USE_ASVS="yes"
USE_SUB_ARGS="yes"
DO_BLASTN_SUB_ARGS="no"

# Info for Kraken analysis
REF_16s_DB="/home/julianna/workspace/gg2_db/mnt/workspace2/jamie/ref/k2__gg2"
THREADS=20

# Input non-file parameters
MAX_SHIFT_PRIMER=2
MAX_MM_PRIMER=2
PRIMER_START_NUM=42
BARCODE_LEN=20
MAX_SHIFT_BARCODE=1
MIN_16S_READS=5
MAX_CONTAM=0.1
FILTER_CORRUPTED= # false
MAX_SHIFT_SUB_ARG=2
MAX_MM_SUB_ARG=0
ALPHA=0.05
BASELINE_GENE="CTX-M"
ARG_THRESHOLD=0.01
MIN_CELLS_PER_ASV=40
FIGURE_DPI=300

# Input MLE parameters
P_MATCH=0.90
P_NONE=0.09
P_ERROR=0.01
ALPHA_PRIOR=1.0
BETA_PRIOR=9.0
MIN_CONFIDENCE=0.95
MIN_NOISE_READS=2
NOISE_CUTOFF_RATIO=0.05

# Intermediate and final output file saving locations

OUT_FOLDER="/home/julianna/workspace/pipeline_output_files_v5/"

FWD_ONLY_16S_FASTQ=$OUT_FOLDER"forward_only_16s.fastq"
REV_ONLY_16S_FASTQ=$OUT_FOLDER"reverse_only_16s.fastq"

KRAKEN_FWD_ONLY_16S_FASTQ=$OUT_FOLDER"kraken_forward_only_16s.fastq"
KRAKEN_REV_ONLY_16S_FASTQ=$OUT_FOLDER"kraken_reverse_only_16s.fastq"

KRAKEN_OUTPUT=$OUT_FOLDER"out.kraken"
KRAKEN_REPORT=$OUT_FOLDER"out.report"

_16S_PACKETS=$OUT_FOLDER"_16s_packets"
ARG_PACKETS=$OUT_FOLDER"arg_packets"
UNCLASSIFIED_PACKETS=$OUT_FOLDER"unclassified_packets"

ARG_REV_FASTQ=$OUT_FOLDER"arg_reverse.fastq"
UNCLASSIFIED_REV_FASTQ=$OUT_FOLDER"unclassified_reverse.fastq"

B_WITH_IDS=$OUT_FOLDER"b_with_ids.txt"
UNFILTERED_BARCODE_SUMMARY_TSV=$OUT_FOLDER"unfiltered_barcode_summary.tsv"
BARCODE_SUMMARY_TSV=$OUT_FOLDER"barcode_summary.tsv"

GLOBAL_ASV_TSV=$OUT_FOLDER"global_asv.tsv"
ASV_BARCODE_SUMMARY_TSV=$OUT_FOLDER"asv_barcode_summary.tsv"

FILTERED_COUNTS_SUMMARY_ARG_TSV=$OUT_FOLDER"filtered_counts_summary_arg.tsv"
FILTERED_BINARY_SUMMARY_ARG_TSV=$OUT_FOLDER"filtered_binary_summary_arg.tsv"
STATS_FILTERING_SUMMARY_ARG_TSV=$OUT_FOLDER"stats_filtering_summary_arg.tsv"

SUB_ARG_SEQS_LIST=$OUT_FOLDER"sub_arg_seqs_list.txt"
SUB_ARG_BARCODE_SUMMARY_TSV=$OUT_FOLDER"sub_arg_barcode_summary.tsv"
EXTRA_MLE_INFO_SUB_ARG_TSV=$OUT_FOLDER"extra_mle_info_sub_arg.tsv"
FILTERED_STATS_CELLS_PER_SUB_ARG_TSV=$OUT_FOLDER"filtered_stats_cells_per_sub_arg.tsv"
FILTERED_SUB_ARG_BARCODE_SUMMARY_TSV=$OUT_FOLDER"filtered_sub_arg_barcode_summary.tsv"

GLOBAL_MLE_TAX_TSV=$OUT_FOLDER"global_mle_tax.tsv"

BLASTN_DB=$OUT_FOLDER"blastn_db"
QUERY_FASTA=$OUT_FOLDER"query_reads.fa"
BLASTN_SUB_ARG_TSV=$OUT_FOLDER"blastn_sub_arg.tsv"

ASV_ARG_TABLE_TSV=$OUT_FOLDER"asv_arg_table.tsv"

# Figures output

ASV_ARG_FIGURE=$OUT_FOLDER"asv_arg_table.png"
BARCODE_GROUP_SIZE_FIGURE=$OUT_FOLDER"barcode_group_size_qc.png"
PRIMER_BALANCE_FIGURE=$OUT_FOLDER"primer_balance_qc.png"


# ---------------------------------------------
# Run the pipeline steps
# ---------------------------------------------

# estimated full pipeline time: 35 mins
echo "Running the full pipeline..."

# ---------------------------------------------

# # Step 0)
# # a few seconds
# echo "[INFO] 0) Validating input files..."
# python validate_inputs.py \
#     --fwd_fastq "$FWD_FASTQ" \
#     --rev_fastq "$REV_FASTQ" \
#     --primers_file "$PRIMER_CSV" \

# # ---------------------------------------------

# # Step 1)
# # <1 min
# echo "[INFO] 1) Extract 16s reads..."
# python extract_16s_reads.py \
#     --fwd_fastq "$FWD_FASTQ" \
#     --rev_fastq "$REV_FASTQ" \
#     --primers_filename "$PRIMER_CSV" \
#     --fwd_only_16s_fastq "$FWD_ONLY_16S_FASTQ" \
#     --rev_only_16s_fastq "$REV_ONLY_16S_FASTQ" \
#     --kraken_fwd_only_16s_fastq "$KRAKEN_FWD_ONLY_16S_FASTQ" \
#     --kraken_rev_only_16s_fastq "$KRAKEN_REV_ONLY_16S_FASTQ" \
#     --max_shift_primer "$MAX_SHIFT_PRIMER" \
#     --max_mm_primer "$MAX_MM_PRIMER" \
#     --primer_start_num "$PRIMER_START_NUM" \

# # ---------------------------------------------

# # Step 2)
# # a few seconds
# echo "[INFO] 2) Run Kraken analysis..."
# kraken2 --use-names --db $REF_16s_DB --threads $THREADS \
#     --report $KRAKEN_REPORT \
#     --paired $KRAKEN_FWD_ONLY_16S_FASTQ $KRAKEN_REV_ONLY_16S_FASTQ \
#     > $KRAKEN_OUTPUT

# # ---------------------------------------------

# # Step 3)
# # 7 mins
# echo "[INFO] 3) Create ID packets..."
# python create_ID_packets.py \
#     --fwd_fastq "$FWD_FASTQ" \
#     --rev_fastq "$REV_FASTQ" \
#     --primers_filename "$PRIMER_CSV" \
#     --kraken_output "$KRAKEN_OUTPUT" \
#     --kraken_report "$KRAKEN_REPORT" \
#     --max_shift_primer "$MAX_SHIFT_PRIMER" \
#     --max_mm_primer "$MAX_MM_PRIMER" \
#     --primer_start_num "$PRIMER_START_NUM" \
#     --barcode_len "$BARCODE_LEN" \
#     --_16s_packet_filename "$_16S_PACKETS" \
#     --arg_packet_filename "$ARG_PACKETS" \
#     --unclassified_packet_filename "$UNCLASSIFIED_PACKETS" \
#     --arg_rev_fastq "$ARG_REV_FASTQ" \
#     --unclassified_rev_fastq "$UNCLASSIFIED_REV_FASTQ" 

# # ---------------------------------------------

# # Step 4)
# # 4 mins
# echo "[INFO] 4) Generate barcode with IDs file..."
# python match_barcodes_to_IDs_revised.py \
#     --_16s_rev_fastq "$REV_ONLY_16S_FASTQ" \
#     --arg_rev_fastq "$ARG_REV_FASTQ" \
#     --unclassified_rev_fastq "$UNCLASSIFIED_REV_FASTQ" \
#     --b_with_ids_filename "$B_WITH_IDS" \
#     --_16s_packet_filename "$_16S_PACKETS" \
#     --arg_packet_filename "$ARG_PACKETS" \
#     --unclassified_packet_filename "$UNCLASSIFIED_PACKETS" \
#     --max_shift_barcode "$MAX_SHIFT_BARCODE" \
#     --barcode_len "$BARCODE_LEN" 

# # ---------------------------------------------

# # Step 5)
# # 22 mins (14 mins making summary + 8 mins filtering)
# echo "[INFO] 5) Generate barcode summary (includes MLE taxonomic classification, ARG counts, and filtering barcodes)..."
# python barcode_summary.py \
#     --b_with_ids_filename "$B_WITH_IDS" \
#     --_16s_packet_filename "$_16S_PACKETS" \
#     --arg_packet_filename "$ARG_PACKETS" \
#     --unfiltered_tsv_filename "$UNFILTERED_BARCODE_SUMMARY_TSV" \
#     --tsv_filename "$BARCODE_SUMMARY_TSV" \
#     --primers_filename "$PRIMER_CSV" \
#     --min_16s_reads "$MIN_16S_READS" \
#     --max_contam "$MAX_CONTAM" \
#     --p_match "$P_MATCH" \
#     --p_none "$P_NONE" \
#     --p_error "$P_ERROR" \
#     --alpha_prior "$ALPHA_PRIOR" \
#     --beta_prior "$BETA_PRIOR" \
#     --min_confidence "$MIN_CONFIDENCE" \
#     --min_noise_reads "$MIN_NOISE_READS" \
#     --noise_cutoff_ratio "$NOISE_CUTOFF_RATIO"

# # ---------------------------------------------

# # Step 6)
# if [ "$USE_ASVS" = "yes" ]; then
#     # Run ASV step of pipeline
#     # a few seconds
#     echo "[INFO] 6) Add ASV information to barcode summary (and filter by ASV status)..."
#     python asv_typing_revised.py \
#         --barcode_summary_tsv "$BARCODE_SUMMARY_TSV" \
#         --b_with_ids "$B_WITH_IDS" \
#         --fwd_16s_fastq "$FWD_ONLY_16S_FASTQ" \
#         --rev_16s_fastq "$REV_ONLY_16S_FASTQ" \
#         --asv_barcode_summary_tsv "$ASV_BARCODE_SUMMARY_TSV" \
#         --global_asv_tsv "$GLOBAL_ASV_TSV" \
#         --primers_file "$PRIMER_CSV" \
#         --filter_corrupted "$FILTER_CORRUPTED"
#     INPUT_FOR_FILTER_ARGS_BARCODE_SUMMARY_TSV=$ASV_BARCODE_SUMMARY_TSV
#     FIRST_GENE_COLUMN_NUM=14
# elif [ "$USE_ASVS" = "no" ]; then
#     echo "[INFO] 6) [SKIP] Skipping ASV step..."
#     INPUT_FOR_FILTER_ARGS_BARCODE_SUMMARY_TSV=$BARCODE_SUMMARY_TSV
#     FIRST_GENE_COLUMN_NUM=5
# fi

# # ---------------------------------------------

# # Step 7)
# # a few seconds
# echo "[INFO] 7) Filter barcode summary by ARGs..."
# python filter_args.py \
#     --input_arg_barcode_summary_tsv "$INPUT_FOR_FILTER_ARGS_BARCODE_SUMMARY_TSV" \
#     --primers_file "$PRIMER_CSV" \
#     --filtered_counts_summary_arg_tsv "$FILTERED_COUNTS_SUMMARY_ARG_TSV" \
#     --filtered_binary_summary_arg_tsv "$FILTERED_BINARY_SUMMARY_ARG_TSV" \
#     --stats_filtering_summary_arg_tsv "$STATS_FILTERING_SUMMARY_ARG_TSV" \
#     --alpha "$ALPHA"

# ---------------------------------------------

# Step 8)
if [ "$USE_SUB_ARGS" = "yes" ]; then
    # Run sub-ARG step of pipeline
    # ~1 min
    echo "[INFO] 8) Create sub-ARG based barcode summary..."
    python sub_arg_database_revised.py \
        --filtered_counts_summary_arg_tsv "$FILTERED_COUNTS_SUMMARY_ARG_TSV" \
        --b_with_ids "$B_WITH_IDS" \
        --arg_packets "$ARG_PACKETS" \
        --fwd_fastq "$FWD_FASTQ" \
        --rev_fastq "$REV_FASTQ" \
        --sub_arg_barcode_summary_tsv "$SUB_ARG_BARCODE_SUMMARY_TSV" \
        --primers_file "$PRIMER_CSV" \
        --sub_arg_seqs_list "$SUB_ARG_SEQS_LIST" \
        --filtered_sub_arg_barcode_summary_tsv "$FILTERED_SUB_ARG_BARCODE_SUMMARY_TSV" \
        --baseline_gene "$BASELINE_GENE" \
        --filtered_stats_cells_per_sub_arg_tsv "$FILTERED_STATS_CELLS_PER_SUB_ARG_TSV" \
        --alpha "$ALPHA" \
        --max_shift_sub_arg "$MAX_SHIFT_SUB_ARG" \
        --max_mm_sub_arg "$MAX_MM_SUB_ARG" \
        --extra_mle_info_sub_arg_tsv "$EXTRA_MLE_INFO_SUB_ARG_TSV" \
        --p_match "$P_MATCH" \
        --p_none "$P_NONE" \
        --p_error "$P_ERROR" \
        --alpha_prior "$ALPHA_PRIOR" \
        --beta_prior "$BETA_PRIOR" \
        --min_confidence "$MIN_CONFIDENCE" \
        --min_noise_reads "$MIN_NOISE_READS" \
        --noise_cutoff_ratio "$NOISE_CUTOFF_RATIO"

    FINAL_BARCODE_SUMMARY_TSV=$FILTERED_SUB_ARG_BARCODE_SUMMARY_TSV

elif [ "$USE_SUB_ARGS" = "no" ]; then
    echo "[INFO] 8) [SKIP] Skipping sub-ARG step..."
    FINAL_BARCODE_SUMMARY_TSV=$FILTERED_COUNTS_SUMMARY_ARG_TSV
fi

# ---------------------------------------------

# Step 9)
if [ "$DO_BLASTN_SUB_ARGS" = "yes" ]; then 
    if [ "$USE_SUB_ARGS" = "yes" ]; then     # must use sub-ARGs in order to have BLASTN sub-ARG table
        echo "[INFO] 9) Create BLASTN sub-ARG table..."
        python blastn_sub_arg.py \
            --sub_arg_seqs_list "$SUB_ARG_SEQS_LIST" \
            --blastn_sub_arg_tsv "$BLASTN_SUB_ARG_TSV" \
            --query_fasta "$QUERY_FASTA" \
            --input_fasta "$DB_SUB_ARG_FASTA" \
            --db "$BLASTN_DB" \
            --final_barcode_summary_tsv "$FINAL_BARCODE_SUMMARY_TSV" \
            --first_gene_column_num "$FIRST_GENE_COLUMN_NUM" 
    fi

elif [ "$DO_BLASTN_SUB_ARGS" = "no" ]; then
    echo "[INFO] 9) [SKIP] Skipping BLASTN sub-ARG table step..."

elif [ "$USE_SUB_ARGS" = "no" ]; then
    echo "[INFO] 9) [SKIP] Skipping BLASTN sub-ARG table step..."

fi

# ---------------------------------------------
















FIRST_GENE_COLUMN_NUM=14























# Step 10)
# a few seconds
echo "[INFO] 10) Create figures (ASV-ARG Table, Barcode Group Size QC, and Primer Balance QC)..."
python figures_program.py \
    --use_asvs_str "$USE_ASVS" \
    --unfiltered_barcode_summary_tsv "$UNFILTERED_BARCODE_SUMMARY_TSV" \
    --final_asv_barcode_summary_tsv "$FILTERED_COUNTS_SUMMARY_ARG_TSV" \
    --asv_barcode_summary_no_sub_args_tsv "$FILTERED_COUNTS_SUMMARY_ARG_TSV" \
    --primers_file "$PRIMER_CSV" \
    --b_with_ids "$B_WITH_IDS" \
    --asv_arg_table_tsv "$ASV_ARG_TABLE_TSV" \
    --global_asv_tsv "$GLOBAL_ASV_TSV" \
    --asv_arg_figure "$ASV_ARG_FIGURE" \
    --barcode_group_size_figure "$BARCODE_GROUP_SIZE_FIGURE" \
    --primer_balance_figure "$PRIMER_BALANCE_FIGURE" \
    --first_gene_column_num "$FIRST_GENE_COLUMN_NUM" \
    --global_mle_tax_tsv "$GLOBAL_MLE_TAX_TSV" \
    --arg_threshold "$ARG_THRESHOLD" \
    --min_cells_per_asv "$MIN_CELLS_PER_ASV" \
    --figure_dpi "$FIGURE_DPI"

# ---------------------------------------------

echo ""
echo ""
echo "[INFO] Completed full pipeline!"
echo ""
echo "    Barcode Summaries:"
echo ""
echo "    - Unfiltered: $UNFILTERED_BARCODE_SUMMARY_TSV"
echo ""
echo "    - Initial Filtering (filter out barcodes with low 16s reads, high contamination,"
echo "       and/or unclassified taxonomy; filter out taxonomic classifications with low # of barcodes):"
echo "       $BARCODE_SUMMARY_TSV"
echo ""
if [ "$USE_ASVS" = "yes" ]; then
    echo "    - ASV Filtering (filter out mixed ASVs, and optionally corrupted single ASVs):" 
    echo "       $ASV_BARCODE_SUMMARY_TSV"
    echo ""
fi
echo "    - ARG Filtering (nullify low ARG read values to filter out ARG noise; based on Poisson distribution):"
echo "       $FILTERED_COUNTS_SUMMARY_ARG_TSV"
echo ""
if [ "$USE_SUB_ARGS" = "yes" ]; then
    echo "    - Sub-ARG Filtering (filter out sub-ARGs with low read values, to preserve only top few sub-ARGs):"
    echo "       $FILTERED_SUB_ARG_BARCODE_SUMMARY_TSV"
    echo ""
fi
echo "    FINAL BARCODE SUMMARY: $FINAL_BARCODE_SUMMARY_TSV"
echo ""
echo ""
echo "    Other Helpful Files:"
echo "    - MLE taxonomic classifications (most to least common):"
echo "      $GLOBAL_MLE_TAX_TSV"
if [ "$USE_ASVS" = "yes" ]; then
    echo "    - ASV names to sequences (most to least common):"
    echo "      $GLOBAL_ASV_TSV"
fi
if [ "$USE_SUB_ARGS" = "yes" ]; then
    echo "    - Sub-ARG names to sequences:"
    echo "      $SUB_ARG_SEQS_LIST"
fi
if [ "$DO_BLASTN_SUB_ARGS" = "yes" ]; then 
    if [ "$USE_SUB_ARGS" = "yes" ]; then
        echo "    - BLASTN sub-ARG table:"
        echo "      $BLASTN_SUB_ARG_TSV"
    fi
fi
echo ""
echo ""
echo "    Figures:"
echo ""
echo "    1) ASV-ARG Table:"
echo ""
echo "    - Numeric Table: $ASV_ARG_TABLE_TSV"
echo "    - Heat Map Figure: $ASV_ARG_FIGURE"
echo ""
echo "    2) Barcode Group Size QC: $BARCODE_GROUP_SIZE_FIGURE"
echo ""
echo "    3) Primer Balance QC: $PRIMER_BALANCE_FIGURE"
echo ""