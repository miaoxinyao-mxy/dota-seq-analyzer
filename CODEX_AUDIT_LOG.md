# CODEX Audit Log

Audit date: 2026-08-28. Repository commit audited: 5b91860. The handoff note is treated as hypotheses, not ground truth. This session implemented only the four requested low-risk engineering fixes; PI and scientific items were not changed.

## 1. Audit summary

| ID | Handoff item | Status | Actual issue | Recommended action | Priority |
| -- | ------------ | ------ | ------------ | ------------------ | -------- |
| 1 | Generalize long target-name wrapping | READY TO FIX | Heat-map code wraps only one literal OXA label | Add generic display-only wrapping | Medium |
| 2 | Generalize Klebsiella correction | PI / SCIENTIFIC DECISION REQUIRED | A hard-coded Klebsiella rewrite exists; broader correctness is unknown | Approve an authoritative normalization policy | High |
| 3 | Validate optional reference FASTA | READY TO FIX | CLI checks only that the path is a file | Define and add minimum FASTA validation | High |
| 4 | Use non-zero exit for missing paths | CONFIRMED | Standalone stages print an error and return 0 | Use sys.exit(1) or equivalent | High |
| 5 | Translate comments | DOCUMENTATION ONLY | Non-English comments remain but do not affect execution | Separate cleanup if desired | Low |
| 6 | TEM_A/TEM_B and separate R1/R2 shifts | ALREADY FIXED | One shift is applied to the concatenated R1|R2 core | Compare R1 and R2 independently, then cluster only when both pass | High |
| 7 | Standardize terminology | DOCUMENTATION ONLY | Internal ARG/sub-ARG terms remain; public CLI is generic | Standardize approved user-facing terms | Medium |
| 8 | Make taxonomy minimum-cell filter optional | READY TO FIX | Stage 2 always runs with min_barcodes=10 | Add an explicit option after choosing the default | Medium |
| 9 | Family vs single target semantics | PI / SCIENTIFIC DECISION REQUIRED | Current Mode supports blank/ssr; legacy single/family is compatibility only | Define biological semantics | High |
| 10 | Zero-inflation in Poisson lambda | PI / SCIENTIFIC DECISION REQUIRED | Global-rate Poisson is used; no zero-inflated model exists | Analyze real data before changing model | High |
| 11 | Add within_cluster_boundary | ALREADY FIXED | Function is present and called in ASV and sub-ARG clustering | No change | — |
| 12 | Update get_genes_eligible_for_sub_args header | OUTDATED | Implementation already uses primer Mode | No code change; wording review only | Low |
| 13 | Use packet JSONL instead of R2 FASTQs | CANNOT VERIFY | Current pipeline intentionally consumes derived R2 FASTQs | Benchmark and trace fields before considering change | Low |
| 14 | Put json.dumps in format_packet | DOCUMENTATION ONLY | Manual JSON construction is fragile but no failure was demonstrated | Optional focused cleanup with serialization tests | Low |
| 15 | qnrB family/single and missing BLAST hits | CANNOT VERIFY | Exact primer/reference fixture is absent from repo evidence | Reproduce with the original files | High |
| 16 | Sub-ARG coordinates 110 vs ASV 120 | PI / SCIENTIFIC DECISION REQUIRED | The two modules use different extraction ends | Confirm intended biological coordinates | High |
| 17 | Shared alpha for Poisson and decay | PI / SCIENTIFIC DECISION REQUIRED | Both default to 0.05 but are separate parameters | Decide whether one or two thresholds are intended | High |
| 18 | Rename CoreASV to ASV | READY TO FIX | Output still emits CoreASV_* | Approve compatibility impact, then rename consistently | Medium |
| 19 | Annotate heat-map values above 1% | READY TO FIX | Heat map has no numeric cell annotations | Add only after display requirements are approved | Medium |
| 20 | Remove Kraken taxid suffixes | PI / SCIENTIFIC DECISION REQUIRED | Labels retain suffixes such as _737866 | Decide whether suffixes are provenance | Medium |
| 21 | Show species for each ASV | PI / SCIENTIFIC DECISION REQUIRED | ASV uses the most specific existing MLE taxonomy; no independent species call | Define evidence and confidence rules | High |
| 23 | ASV also shares one shift across paired cores | ALREADY FIXED | ASV comparison uses one shift on concatenated R1/R2 | Add pair-aware independent-shift comparison | High |
| 22 | GreenGenes2/KrakenTools citations | DOCUMENTATION ONLY | Paper attribution, not pipeline correctness | Confirm citations in manuscript materials | Medium |

## 2. Detailed investigation notes

### Item 1 — Heat-map target-name formatting

**Handoff claim**

Only OXA-48 like (oxa-3) is split; long names should be handled generally.

**Files/functions inspected**

src/dota_seq_analyzer/figures_program.py, make_asv_arg_table().

**What the code actually does / evidence**

Lines 95–97 compare each label to the single literal OXA-48 like (oxa-3) and replace its space with a newline. Other labels are unchanged.

**Verdict**

Confirmed narrow presentation limitation; not a data-assignment bug.

**Recommended action / risk**

Generic display-only wrapping is reasonable, but line-break width and punctuation policy are unverified. Preserve TSV/JSON target values.

**Assumptions / unverified points**

Desired figure width and wrapping style are unspecified.

**Changes made**

None — audit only.

### Item 2 — Klebsiella taxonomic correction

**Handoff claim**

Klebsiella is corrected and other taxa may need correction.

**Files/functions inspected**

src/dota_seq_analyzer/filter_barcodes.py, filter_barcodes_in_df().

**What the code actually does / evidence**

Surviving rows containing G - Klebsiella are rewritten to Klebsiella pneumoniae complex with species None. This is a literal special case.

**Verdict**

The special case is confirmed. Broader or different corrections cannot be inferred from code.

**Recommended action / risk**

PI / SCIENTIFIC DECISION REQUIRED. A normalization table could alter reported taxonomy and must not be guessed.

**Assumptions / unverified points**

Authoritative database naming and desired corrections were not supplied.

**Changes made**

None — audit only.

### Item 3 — Reference FASTA validation

**Handoff claim**

Optional reference FASTA should be validated by input validation.

**Files/functions inspected**

src/dota_seq_analyzer/cli.py, validate_inputs.py, blastn_sub_arg.py.

**What the code actually does / evidence**

cli.py checks reference.is_file() only. The reference is not passed to validate_inputs and no FASTA record/sequence validation occurs before makeblastdb.

**Verdict**

Confirmed validation gap.

**Recommended action / risk**

READY TO FIX after defining the minimum contract: readable FASTA, at least one record, non-empty IDs/sequences, and an approved alphabet. Incorrect restrictions could reject valid references.

**Assumptions / unverified points**

Ambiguous bases, duplicate IDs, and non-nucleotide characters are unresolved.

**Changes made**

None — audit only.

### Item 4 — Missing-path exit status

**Handoff claim**

Missing paths should stop with sys.exit(1), not return.

**Files/functions inspected**

main() in extract_16s_reads.py, create_ID_packets.py, match_barcodes_to_IDs_revised.py, barcode_summary.py, asv_typing_revised.py, filter_args.py, sub_arg_database_revised.py, and blastn_sub_arg.py.

**What the code actually does / evidence**

These modules use print(error) followed by bare return in missing-input branches. The public CLI uses subprocess.run(check=True), so a child returning 0 can permit invalid downstream state.

**Verdict**

Confirmed bug.

**Recommended action / risk**

Change only input-error branches to non-zero termination and add a subprocess failure test. Preserve successful returns and output behavior.

**Assumptions / unverified points**

The exact mechanism may be sys.exit(1) or parser.error; non-zero termination is the required behavior.

**Changes made**

None — audit only.

### Item 5 — Comments and terminology

**Handoff claim**

Comments should be English and terminology standardized.

**Files/functions inspected**

mle_revised.py, sub_arg_database_revised.py, CLI/help text and output labels.

**What the code actually does / evidence**

Some comments contain Chinese text/emoji. Internal names and columns still contain ARG/sub-ARG, while the public entry point describes generic targeted genes and phase variation.

**Verdict**

Documentation/readability only; no runtime bug demonstrated.

**Recommended action / risk**

Do separately after terminology is approved. Renaming columns can break downstream users.

**Changes made**

None — audit only.

### Items 6–22 — Evidence and verdicts

**Files/functions inspected**

sub_arg_database_revised.py, helper_functions.py, filter_args.py, filter_barcodes.py, asv_typing_revised.py, match_barcodes_to_IDs_revised.py, create_ID_packets.py, figures_program.py, blastn_sub_arg.py, cli.py, export_results.py.

**What the code actually does / evidence**

* Item 6: sub-ARG cores are stored as one R1|R2 string and passed to one semi_global_distance() call, so one shift is shared by both reads. This is confirmed as an implementation bug under the clarified design.
* Item 7: public naming is generic but internal terminology remains.
* Item 8: filter_barcodes_in_df always applies Stage 2 with min_barcodes=10; CLI has no switch.
* Item 9/12: get_target_modes reads blank or ssr; legacy single/family is accepted only for legacy validation. get_genes_eligible_for_sub_args uses Mode, or all targets for reference reconstruction.
* Item 10: filter_args calculates r_j from global target reads/global 16S reads and lambda_ij from cell 16S reads; it is an ordinary Poisson calculation.
* Item 11: within_cluster_boundary exists and is called in ASV and sub-ARG clustering.
* Item 13: match_barcodes_to_IDs_revised.py consumes derived R2 FASTQs; packet fields and replacement performance were not benchmarked.
* Item 14: format_packet manually constructs JSON-like text; no failing fixture was demonstrated.
* Item 15: exact qnrB inputs are not present in the inspected evidence.
* Item 16: ASV coordinates are [30:120] and [70:120]; sub-ARG coordinates are [30:110] and [70:110].
* Item 17: filter_args and sub-ARG defaults are both 0.05, but separate parameters; top-level CLI does not expose them.
* Item 18: ASV table/header uses CoreASV_* and Core_ASV_ID.
* Item 19: heat-map values are thresholded/plotted without text annotations.
* Item 20: taxonomy strings retain suffixes such as Enterobacterales_737866.
* Item 21: ASV labels use the most specific MLE taxonomy available; no independent species assignment is performed.

**Verdict**

Items 9, 10, 16, 17, 20, and 21 require PI/scientific decisions. Item 6 is a confirmed implementation bug under the clarified intended design. Items 7, 13, and 14 are not demonstrated bugs. Items 8, 18, and 19 are real feature/behavior gaps. Item 12 is already reflected in implementation. Item 15 cannot be reproduced from this repository.

**Recommended action / risk**

Do not change scientific behavior until the relevant decision is recorded. For implementation gaps, add the smallest targeted change and validation only after the user approves the behavior.

**Assumptions / unverified points**

No original primer/reference fixture, formal taxonomy normalization table, statistical model specification, or publication display specification was supplied.

**Changes made**

None — audit only.

## 3. Confirmed bugs

1. Standalone missing-input paths can return status 0 after printing an error.
2. Heat-map wrapping is hard-coded to one target label.
3. Optional reference FASTA is checked for existence but not content.
4. ASV and sub-ARG paired-core comparison applies one shift to both R1 and R2, inconsistent with the clarified intended design.

## 4. Incorrect or outdated handoff notes

* within_cluster_boundary is already implemented and used.
* Current sub-ARG eligibility uses primer Mode; the single/family claim is outdated for the current implementation.
* The requested get_genes_eligible_for_sub_args behavior is already implemented.
* Public R1/R2 terminology is already used, with internal compatibility aliases retained.

## 5. PI / scientific decisions

Taxonomic normalization; family vs single semantics; zero-inflated modeling; extraction coordinates; shared versus separate alpha; taxid display/provenance; and species evidence for ASVs. The TEM/ASV independent-shift issue is confirmed implementation behavior, not a PI decision.

## 6. Code changes made

* validate_inputs.py: added check_reference_fasta() for structural FASTA validation; missing-input main() branches now exit non-zero.
* cli.py: calls check_reference_fasta() when -r/--reference is supplied.
* extract_16s_reads.py, create_ID_packets.py, match_barcodes_to_IDs_revised.py, barcode_summary.py, asv_typing_revised.py, filter_args.py, sub_arg_database_revised.py, blastn_sub_arg.py, and figures_program.py: missing required path branches now terminate non-zero; successful paths are unchanged.
* figures_program.py::make_asv_arg_table(): replaced the single OXA display exception with generic wrapping and added numeric text to retained non-zero heat-map cells. Original target values and matrix thresholding remain unchanged.

**Behavior before**

Missing-input child scripts could return status 0; reference files were checked only for existence; one target label had a hard-coded line break; heat-map cells had color only.

**Behavior after**

Missing-input child scripts return status 1; supplied references must contain structurally valid FASTA records; all target labels are display-wrapped at width 18; retained heat-map values are printed to two decimals.

**Remaining assumptions**

FASTA validation is intentionally structural and does not impose a nucleotide alphabet, uniqueness policy, or gene/family interpretation. Heat-map width 18 and two-decimal display are presentation choices.

## 7. Tests / validation performed

* python -m py_compile src/dota_seq_analyzer/*.py: PASS.
* Focused check_reference_fasta tests for valid, empty-record, sequence-before-header, and blank-header cases: PASS.
* Focused filter_args.py missing-input subprocess test: exit code 1 and expected error: PASS.
* Package CLI --help and invalid-reference rejection: PASS; a package-import compatibility issue was found and fixed within the reference-validation change.
* Figures stage run in /home/xy/02_dota-seq/test_runs/audit_figures with successful test05 inputs: PASS; all three PNGs generated.
* CoreASV naming, Stage 2 filtering, and PI/scientific items were not changed.

## 8. Open items / next actions

1. Review and commit/push the completed four scoped engineering fixes.
2. Decide whether structural FASTA validation should later add alphabet or duplicate-ID rules.
3. Approve heat-map label width and numeric precision if publication formatting requires it.
4. Obtain PI decisions for the scientific items listed above.
5. Only then consider ASV renaming and packet-input optimization.

## 9. Session history

### 2026-08-28

Inspected and re-inspected the current repository against the undergraduate handoff. Implemented and focused-tested only the four requested low-risk engineering fixes. CoreASV naming and Stage 2 filtering were traced but not changed; PI/scientific items remain untouched. Next session should review the diff and decide whether to commit/push these scoped changes.

### Item 23 — ASV independent-shift comparison — 2026-08-28

**Handoff clarification**

R1 and R2 are intended to allow independent valid shifts; TEM_A/TEM_B is only a test case exposing the general rule.

**Files/functions inspected**

src/dota_seq_analyzer/asv_typing_revised.py: extract_core(), semi_global_distance(), within_cluster_boundary(), summarize_barcode(). Downstream uses: global_asv.tsv generation, Assigned_core_asv, figures_program.py, and export_results.py.

**Current behavior**

extract_core() returns one R1|R2 string. summarize_barcode() and within_cluster_boundary() each call semi_global_distance() once on that concatenated string, so a single shift is applied to both components.

**Evidence**

A synthetic pair with independent valid shifts gave distance 0 for R1, distance 0 for R2, and distance 7 when concatenated and passed through the current ASV distance function.

**Root cause**

This is an accidental consequence of concatenating paired reads before applying a single-read distance function, not an explicit biological rule.

**Verdict**

Confirmed implementation bug under the clarified original design. It affects ASV clustering generally, not only TEM_A/TEM_B.

**Recommended action**

Do not modify until separately authorized. The minimum fix is a pair-aware distance helper that compares the R1 and R2 components independently, then uses it in candidate and boundary comparisons. Preserve ASV thresholds, statuses, identifiers, and output schema.

**Expected downstream effects**

Some reads currently assigned to different ASVs may cluster together. This can change per-cell ASV status, global ASV membership, heat-map grouping, and JSONL ASV assignments. Blank-Mode target detection is unaffected.

**Assumptions / unverified points**

The synthetic check proves the shared-shift implementation behavior, not the final mismatch aggregation policy. The original TEM fixture should be used for regression validation.

**Changes made**

None — audit only.

## Session history — follow-up

### 2026-08-28

Re-audited ASV after the original developer clarified that R1 and R2 must allow independent valid shifts. Confirmed the same shared-shift implementation bug as sub-ARG comparison. A synthetic read pair produced component distances 0 and 0 but concatenated distance 7. No ASV code was changed. Next action is to design and separately authorize the minimum pair-aware distance fix, with the original TEM fixture as regression validation.

## Resolution update — independent R1/R2 shifts — 2026-08-28

**Files/functions changed**

* src/dota_seq_analyzer/asv_typing_revised.py: added paired_semi_global_distance(); summarize_barcode() and within_cluster_boundary() now use it for paired cores.
* src/dota_seq_analyzer/sub_arg_database_revised.py: added paired_semi_global_distance(); cluster_sub_arg_seqs() and within_cluster_boundary() now use it.

**Behavior before**

Each module stored paired cores as R1|R2 and applied one semi_global_distance() shift to the concatenated string, forcing R1 and R2 to share a shift.

**Behavior after**

R1 and R2 are split and compared independently, each allowing the existing max shift. The pair distance is the larger component distance, so both components must satisfy the existing mismatch threshold. Extraction coordinates, thresholds, cluster ordering, ASV identifiers, statuses, and output formats remain unchanged.

**Validation**

python -m py_compile src/dota_seq_analyzer/*.py: PASS. Synthetic pair with R1 distance 0 and R2 distance 0 under independent shifts: PASS for both helpers. ASV summarize_barcode() clustered the pair: PASS. sub-ARG cluster_sub_arg_seqs() clustered the pair: PASS.

**Remaining assumptions**

The mismatch threshold is treated as applying independently to each component; this matches the clarified design and existing zero-mismatch default. Original TEM fixture and a real-data end-to-end rerun remain recommended. No PI/scientific items were changed.


## Follow-up validation — TEM_A/TEM_B-like regression — 2026-08-28

The original TEM_A/TEM_B primer/reference fixture was not present in the repository or existing test directories. A fixture with the same independent-shift pattern was therefore passed through the actual `sub_arg_database_revised.py::cluster_sub_arg_seqs()` implementation.

**Result**

* Old shared-shift distance: 7; sequences would remain separate.
* New independent-shift distance: 0.
* Resulting cluster count: 1.
* Regression status: PASS.

This validates the implementation path for the reported behavior, but does not replace validation against the original biological TEM_A/TEM_B sequences.
