# DoTA-Seq Analyzer canonical result schema v3.0.0

This document is the frozen contract for the first canonical result set. Every
pipeline invocation generates one UUIDv4 `run_id`; the same value is required in
all canonical records from that invocation.

JSONL files contain one JSON object per line. Fields not listed as optional are
required. Counts are non-negative integers unless stated otherwise. Numeric QC
values are finite JSON numbers or `null` where the existing pipeline has no
value.

## `cells.jsonl`

One record per final surviving cell/barcode.

```json
{
  "schema_version": "3.0.0",
  "record_type": "cell",
  "run_id": "UUIDv4",
  "cell_barcode": "string...",
  "taxonomy": {
    "ranks": {
      "domain": "Bacteria",
      "phylum": "Bacillota_I",
      "class": "Bacilli_A",
      "order": null,
      "family": null,
      "genus": null,
      "species": null
    },
    "confidence": 1.0,
    "contamination": 0.01
  },
  "taxonomy_evidence": {
    "total_16s_reads": 12,
    "technical_noise_reads": 1
  },
  "asv": {
    "asv_id": "ASV_1",
    "status": "single_ASV",
    "reads_used": 11,
    "raw_unique_core_sequences": 2,
    "dominant_raw_read_count": 10,
    "coexisting_2bp_reads": 1,
    "unauthorized_secondary_reads": 0,
    "final_cell_asv_reads": 11,
    "max_internal_distance": 1
  },
  "targets": [
    {
      "target_name": "TEM",
      "raw_read_count": 4,
      "filtered_read_count": 3,
      "assignment_type": "sequence_cluster",
      "sequence_cluster_id": "TEM_seq_1",
      "sequence_call_qc": {
        "confidence": 0.99,
        "contamination": 0.01,
        "total_reads": 3,
        "technical_noise_reads": 0
      }
    }
  ]
}
```

Required fields are all fields shown above except:

- `sequence_cluster_id`: required only when `assignment_type` is
  `sequence_cluster`; forbidden otherwise.
- `sequence_call_qc`: optional and present only when the existing reconstructed
  target-sequence MLE output contains QC values for this cell and target.

`targets` is sparse: a target is included only when its raw read count is
greater than zero. `assignment_type` has exactly these meanings:

- `filtered_out`: raw evidence existed, but the filtered count is zero.
- `detected`: filtered evidence remains and no sequence reconstruction was
  assigned or requested.
- `sequence_unresolved`: filtered evidence remains, but reconstruction retained
  only the target parent assignment.
- `sequence_cluster`: filtered evidence was assigned to a reconstructed
  sequence cluster.

`asv.asv_id` is a foreign key to `asvs.jsonl`.
`targets[].sequence_cluster_id` is a foreign key to
`target_sequences.jsonl`.

## `asvs.jsonl`

One record per ASV referenced by a final surviving cell.

```json
{
  "schema_version": "3.0.0",
  "record_type": "asv",
  "run_id": "UUIDv4",
  "asv_id": "ASV_1",
  "core_sequence": {
    "r1": "ACGT...",
    "r2": "TGCA..."
  },
  "final_surviving_cell_count": 42
}
```

All fields are required. `asv_id` follows `ASV_<positive integer>`.
`final_surviving_cell_count` is recalculated from references in `cells.jsonl`.
Canonical v3.0.0 deliberately does not store dominant taxonomy in ASV records.

## `target_sequences.jsonl`

Generated only when target sequence reconstruction is performed. It contains
one record per reconstructed sequence cluster referenced by a final surviving
cell. The file may be empty when reconstruction ran but no sequence cluster
survived.

```json
{
  "schema_version": "3.0.0",
  "record_type": "target_sequence",
  "run_id": "UUIDv4",
  "sequence_cluster_id": "TEM_seq_1",
  "parent_target": "TEM",
  "core_sequence": {
    "r1": "ACGT...",
    "r2": "TGCA..."
  },
  "final_surviving_cell_count": 7,
  "phase_variation": {
    "mode": "ssr",
    "call": "...",
    "motif": "...",
    "repeat_count": 4,
    "start": 10,
    "end": 21,
    "direct_similarity": 0.95,
    "evidence": "..."
  },
  "reference_matches": [
    {
      "subject": "reference identifier",
      "r1_start": 1,
      "r1_end": 90,
      "r1_percent_identity": 100.0,
      "r2_start": 91,
      "r2_end": 140,
      "r2_percent_identity": 100.0
    }
  ]
}
```

All top-level fields are required except:

- `phase_variation`: optional; copied from existing SSR phase-variation output
  when a call exists for the sequence cluster.
- `reference_matches`: optional; present when reference matching was requested,
  and may be an empty array.

The values inside `phase_variation` are the existing output fields; nullable
values remain `null`. `sequence_cluster_id` follows
`<parent_target>_seq_<positive integer>`. `parent_target` must occur in the run's
target panel. `final_surviving_cell_count` is recalculated from final
`cells.jsonl` references and never trusted from an earlier sequence-list count.

## `run_summary.json`

Exactly one JSON object per run.

Required top-level fields:

- `schema_version`: exactly `3.0.0`.
- `record_type`: exactly `run_summary`.
- `run_id`: the run UUIDv4.
- `status`: exactly `completed` for a promoted canonical result set.
- `software`: `name`, `version`, and nullable `git_commit`.
- `inputs`: resolved paths for `r1`, `r2`, `primers`, `taxonomy_database`, and
  nullable `reference_fasta`.
- `target_panel`: ordered records containing `target_name` and effective `mode`.
- `parameters`: the effective runtime configuration described below.
- `read_qc`: `raw_read_pairs`, `passed_read_pairs`, `filtered_read_pairs`, and
  mutually exclusive `rejected` counts from the read-QC stage.
- `read_classification`: actual counts of `accepted_16s_reads`, `target_reads`,
  and `unclassified_reads` in the mutually exclusive packet files.
- `16s_r1_primer_starts`: accepted 16S read counts keyed by each configured
  valid start.
- `barcode_funnel`: direct stage counts for `raw_unique_barcodes`,
  `clustered_barcodes`, `after_stage1_taxonomy_filter`,
  `after_minimum_cells_per_taxon`, `before_asv_filter`,
  `after_asv_status_filter`, `asv_taxonomy_conflicts_removed`, and
  `final_cells`.
- `asv_summary`: `final_asv_count`.
- `target_filtering_summary`: one record per target containing `target_name`,
  `original_positive_cells`, `filtered_out_cells`,
  `remaining_positive_cells`, and `retention_percent`.

Optional top-level field:

- `phase_variation_summary`: present when SSR phase-variation analysis ran;
  contains `cell_calls` and counts keyed by existing call value.

`parameters` records only effective values obtained from actual runtime
arguments or imported from the same constants used by the algorithm. Its
required sections are:

- `analysis_workers`
- `read_qc`: `minimum_read_length`, `minimum_mean_phred`, `barcode_length`,
  `minimum_barcode_q25_bases`
- `primer_matching`: `maximum_shift`, `maximum_mismatches`,
  `r2_primer_start`, `valid_16s_r1_starts`
- `barcode_clustering`: `barcode_length`, `maximum_shift`
- `cell_taxonomy_filtering`: `minimum_16s_reads`,
  `maximum_contamination`, `minimum_cells_per_taxon`
- `taxonomy_mle`: `p_match`, `p_none`, `p_error`, `alpha_prior`, `beta_prior`,
  `minimum_confidence`, `minimum_noise_reads`, `noise_cutoff_ratio`
- `asv`: `r1_start`, `r1_end`, `r2_start`, `r2_end`, `maximum_distance`,
  `maximum_shift`, `minimum_reads`, `mixed_ratio_threshold`,
  `filter_corrupted_single_asv`, `taxonomy_conflict_minimum_cells`,
  `taxonomy_conflict_dominant_phylum_fraction`
- `target_background_filtering`: `alpha`
- `target_sequence_reconstruction`: `performed`, `maximum_shift`,
  `maximum_mismatches`, `alpha`, `r1_start`, `r1_end`, `r2_start`, `r2_end`,
  `include_all_targets`, and an `mle` object with the same fields as
  `taxonomy_mle`

## Mapping from existing final scientific fields

| Existing source | Canonical destination |
| --- | --- |
| `cell_target_matrix.tsv`: Barcode | `cell_barcode` |
| Predicted taxonomy | `taxonomy.ranks` (parsed explicit rank fields) |
| Confidence / Contamination | `taxonomy.confidence` / `taxonomy.contamination` |
| Total # of 16s reads / Technical noise count | `taxonomy_evidence` |
| Assigned_core_asv / Status | `asv.asv_id` / `asv.status` |
| Existing per-cell ASV QC columns | corresponding fields under `asv` |
| `asv_barcode_summary.tsv`: raw target counts | `targets[].raw_read_count` |
| `filtered_counts_summary_arg.tsv`: filtered target counts | `targets[].filtered_read_count` |
| Final target value | `assignment_type` and optional `sequence_cluster_id` |
| Existing reconstructed-sequence MLE table | optional `sequence_call_qc` |
| `global_asv.tsv`: paired core sequence | `asvs.jsonl.core_sequence` |
| `sub_arg_seqs_list.txt`: parent and paired sequence | `target_sequences.jsonl` IDs, parent, and core sequence |
| Existing SSR output | optional `phase_variation` |
| Existing BLAST/reference output | optional `reference_matches` |

No report or intermediate file becomes canonical merely because it is used as
an input to this Phase 1 exporter.

## Validation and publication invariants

Before publication of the result set:

1. Every canonical record has schema version `3.0.0` and the same canonical
   UUIDv4 `run_id`.
2. Cell barcodes, ASV IDs, and target sequence IDs are unique in their files.
3. Every cell ASV foreign key resolves, and no unreferenced ASV record exists.
4. Every cell target-sequence foreign key resolves, and no unreferenced target
   sequence record exists.
5. Each ASV and target-sequence `final_surviving_cell_count` exactly equals its
   number of references in `cells.jsonl`.
6. `passed_read_pairs` equals `accepted_16s_reads + target_reads +
   unclassified_reads`.
7. The sum of 16S primer-start counts equals `accepted_16s_reads`.
8. `barcode_funnel.final_cells` equals the number of cell records, and
   `asv_summary.final_asv_count` equals the number of ASV records.

Canonical files are first written under a run-specific temporary path, reloaded,
and validated. They are promoted only after all schema and relationship checks
pass. The complete run directory is then atomically renamed to the requested
output path, so a failed run cannot leave a partially finalized canonical result
set at that path.
