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
conda activate dota-seq-analyzer
python -m pip install -e .
```

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
dota-seq-analyzer -1 reads1.fastq -2 reads2.fastq -p primers.csv -o results
```

To annotate reconstructed target sequences with a reference FASTA:

```bash
dota-seq-analyzer -1 reads1.fastq -2 reads2.fastq -p primers.csv -r reference.fa -o results
```

An optional AMR reference is included in `database/amr-reference-db.tar.gz`.

Run `dota-seq-analyzer --help` for database overrides, the Kraken2 thread option, and filtering options. By default, taxa must be represented by at least 10 cells. Use `--min-cells-per-taxon 5` to retain taxa with at least five cells, `--min-cells-per-taxon 1` to retain every taxon represented after Stage 1, or `--min-cells-per-taxon 0` to disable Stage 2 taxonomy count filtering.

## Output

The primary results are `dota_seq_analyzer_results.jsonl` and `reports/cell_target_matrix.tsv`. Phase-variation calls are written to `reports/cell_phase_variation.tsv`, and optional BLAST matches to `reports/reference_matches.tsv`.

Intermediate files are written to `tmp/`, report tables to `reports/`, and figures to `figures/`.

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
