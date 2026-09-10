# DoTA-Seq Analyzer

*Single-cell profiling of targeted microbial genes and phase variation*

DoTA-Seq Analyzer links primer-defined targets to individual bacterial cells and their taxonomic classifications in single-cell DoTA-Seq data.

## Installation

```bash
git clone https://github.com/miaoxinyao-mxy/dota-seq-analyzer.git
cd dota-seq-analyzer
git lfs install
git lfs pull
conda env create -f environment.yml
tar -xzf database/dota-seq-analyzer-taxonomy-db.tar.gz -C database
conda run -n dota-seq-analyzer python -m pip install -e .
conda run -n dota-seq-analyzer dota-seq-analyzer --help
```

The bundled GreenGenes2-based Kraken2 taxonomy database is tracked with Git LFS,
extracted by the command above, and discovered automatically by the pipeline.
Use `--taxonomy-db` only to select a different extracted Kraken2 database.

If your interactive shell has already been initialized for Conda, you may
optionally run `conda activate dota-seq-analyzer` before using the software.

## Primer file

The primer CSV contains four columns:

```csv
Primer,F,R,Mode
16s,F_PRIMER,R_PRIMER,
target_1,F_PRIMER,R_PRIMER,
target_2,F_PRIMER,R_PRIMER,ssr
target_3,F_PRIMER,R_PRIMER,ssr
```

Leave `Mode` blank for standard target detection. Use `ssr` for phase-variation targets.

## Run

```bash
conda run -n dota-seq-analyzer dota-seq-analyzer -1 reads1.fastq -2 reads2.fastq -p primers.csv -o results -@ 8
```

Analysis is serial by default. For larger datasets, use `-@ N` or `--threads N` to parallelize the sequential analysis stages across CPU cores; 8 is a reasonable starting point on a multi-core workstation.

After running `conda activate dota-seq-analyzer`, you may omit the
`conda run -n dota-seq-analyzer` prefix from pipeline commands.

To annotate reconstructed target sequences with a reference FASTA:

```bash
conda run -n dota-seq-analyzer dota-seq-analyzer -1 reads1.fastq -2 reads2.fastq -p primers.csv -r reference.fa -o results
```

An optional AMR reference is included in `database/amr-reference-db.tar.gz`.

Run `dota-seq-analyzer --help` for database overrides, parallel analysis, and filtering options. By default, taxa must be represented by at least 10 cells. Use `--min-cells-per-taxon 5` to retain taxa with at least five cells, `--min-cells-per-taxon 1` to retain every taxon represented after Stage 1, or `--min-cells-per-taxon 0` to disable Stage 2 taxonomy count filtering.

## Output

The authoritative v3.0.0 results are:

- `cells.jsonl`: one canonical record per surviving cell, including taxonomy, ASV assignment, QC fields, and target calls.
- `asvs.jsonl`: canonical ASV identifiers, core sequences, and final surviving-cell counts.
- `run_summary.json`: run identity, inputs, parameters, read/cell filtering counts, and result totals.
- `target_sequences.jsonl`: canonical reconstructed target-sequence clusters when sequence reconstruction is performed; otherwise this file is not generated.

These files share one `run_id` and are schema- and relationship-validated before the completed output directory is published. The schema contract is documented in [`docs/canonical-schema-v3.0.0.md`](docs/canonical-schema-v3.0.0.md).

`dota_seq_analyzer_results.jsonl` is retained as a legacy compatibility export. Files under `reports/`, including `cell_target_matrix.tsv`, are derived human-readable tables rather than the canonical source of truth. Phase-variation calls are written to `reports/cell_phase_variation.tsv`, and optional BLAST matches to `reports/reference_matches.tsv`.

Human-readable derived tables are written to `reports/`, and plots are written to `figures/`. Computational intermediates are removed after a successful run by default. Use `--keep-tmp` to preserve the complete `tmp/` directory for development or debugging; failed runs retain their incomplete staging directory and intermediates.

## Optional: PrimerPicker

```bash
python PrimerPicker/primer_picker.py targets.fa --outdir primer_picker_results --seed 123
```

The main output is `primer_picker_results/top-primer-sets.tsv`.

## Citation

DoTA-Seq Analyzer was developed for analysis of single-cell targeted sequencing data based on the DoTA-Seq framework. For the underlying DoTA-Seq method, please cite:

Lan F, Saba J, Ross TD, Zhou Z, Krauska K, Anantharaman K, Landick R, Venturelli OS. Massively parallel single-cell sequencing of diverse microbial populations. *Nature Methods* 21, 228–235 (2024). https://doi.org/10.1038/s41592-023-02157-7

## Contact / Support

If you encounter any bugs or have questions about the code, please feel free to reach out to us:
* **Xinyao Miao**: miaoxinyao.xjtu@gmail.com
* **Julianna Yuen**: juliannayuen07@gmail.com

For general academic inquiries, please contact **Freeman Lan**: freeman.lan@utoronto.ca

## License

This project is released under the MIT License. See [LICENSE](LICENSE).
