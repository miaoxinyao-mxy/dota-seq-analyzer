# Phase 2A report-source audit

Canonical v3 is the authoritative result. This audit records how current TSV
outputs map to canonical data before any legacy report writer is replaced.

| Report content | Current writer | Current upstream sources | Columns / row identity | Ordering and missing values | Canonical recovery |
| --- | --- | --- | --- | --- | --- |
| `reports/cell_target_matrix.tsv` | `sub_arg_database_revised.create_sub_arg_barcode_summary()` when reconstruction runs; otherwise `cli.main()` copies the filtered target table | final ASV cell table, filtered target counts, optional reconstructed target assignments | One row per `Barcode`; 15 fixed taxonomy/ASV fields followed by one column per target in primer-panel order | Cell order is inherited from the filtered cell table. Absent/filtered targets are `0`; unresolved reconstructed targets are `<target>_parent`; resolved calls are `<target>_seq_N` | Lossless from ordered `cells.jsonl` plus target-panel order in `run_summary.json` |
| ASV summary (`tmp/global_asv.tsv`) | `asv_typing_revised.conduct_asv_typing()` | surviving ASV cell assignments and global core sequences | `Core_ASV_ID`, `cell_count`, `core_sequence`; one row per surviving ASV | Numeric ASV-ID order inherited from ASV creation; no missing values | Lossless from `asvs.jsonl` |
| Taxonomy summary (`tmp/global_mle_tax.tsv`) | `figures_program.write_global_tax_classification_file()` | final cell matrix | `MLE_Taxonomic_Classification`, `Cell_Count`; one row per complete lineage string | Count descending; ties retain first-cell occurrence order. Unclassified ranks are literal `None` inside the lineage | Lossless from ordered `cells.jsonl` |
| Target summary (`tmp/stats_filtering_summary_arg.tsv`) | `filter_args.filter_args()` | raw and Poisson-filtered per-cell target counts | `ARG`, `Original`, `Filtered Out`, `Remaining`, `% Retention`; one row per panel target | Primer-panel order; zero-original targets have `0.0` retention | Lossless from `run_summary.json.target_filtering_summary` |
| `reports/cell_phase_variation.tsv` | `phase_variation.write_cell_calls()` | reconstructed sequence list, final cell matrix, primer modes | `Barcode`, `Target`, `Mode`, `Sequence_assignment`, `Call`, `Motif`, `Repeat_count`, `Start`, `End`, `Direct_similarity`, `Evidence`; one row per cell/SSR-target call | Cell order followed by panel-target order. Unavailable optional call values are empty TSV fields | Lossless by joining `cells.jsonl` sequence-cluster foreign keys to `target_sequences.jsonl.phase_variation` |
| `reports/reference_matches.tsv` | `blastn_sub_arg.make_summary_table()` | reconstructed sequence list, final assignments, BLAST results | `sub-ARG_name`, `subject`, R1/R2 coordinates and identities, `read1_seq`, `read2_seq`; one row per sequence-cluster/reference match | Sequence-list order; subjects alphabetic within a sequence. No defined missing values | Lossless from `target_sequences.jsonl.reference_matches` plus its R1/R2 core sequence |

The exact fixed columns of `cell_target_matrix.tsv`, before the dynamic target
columns, are:

```text
Barcode
Predicted taxonomy
Confidence
Contamination
Total # of 16s reads
Technical noise count
Reads_used_for_ASV
Raw unique_core_sequences
Dominant_raw_read_count
Coexisting_2bp_reads
Unauthorized_secondary_reads
Final_cell_asv_reads
Max_internal_distance
Assigned_core_asv
Status
```

No other TSV is written under `reports/` by the current production CLI. TSVs
under `tmp/` other than the three summaries above remain computational or
figure/QC intermediates and are outside Phase 2A migration.

`summary.tsv` is a new convenience report derived only from existing
`run_summary.json` fields. It does not define or calculate new scientific
metrics.
