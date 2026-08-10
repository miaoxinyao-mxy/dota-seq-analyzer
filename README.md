# ARGMapper: AMR mapping for DoTA-Seq

This repository contains ARGMapper, a DoTA-Seq workflow for mapping antibiotic-resistance genes (ARGs) to bacterial taxonomic classifications and visualizing AMR patterns.

PrimerPicker is included as an optional utility for designing multiplex PCR primer pools. It is not required for running ARGMapper.

## Repository layout

```text
src/argmapper/       ARGMapper pipeline modules and driver configuration
PrimerPicker/        Optional primer-design utility and original notebooks
pyproject.toml        Python package metadata
LICENSE               MIT license
```

## ARGMapper

ARGMapper is intended for paired-end, short-read bacterial samples generated as part of the DoTA-Seq workflow. It uses 16S taxonomic classification, barcode-level processing, ASV information, ARG filtering, optional sub-ARG analysis, and statistical classification to produce per-cell summaries and AMR/ASV–ARG maps.

### Installation

Create a Python environment with Python 3.13 or newer and install the package from the repository root:

```bash
python -m pip install .
```

The full pipeline also depends on external DoTA-Seq tools and databases, including Kraken2 and, when sub-ARG analysis is enabled, BLAST. Configure their paths and analysis parameters in `src/argmapper/driver.sh` before running the workflow.

### Required inputs

The standard workflow expects:

- forward and reverse FASTQ files
- a DoTA-Seq primer CSV
- a sub-ARG database FASTA when sub-ARG analysis is enabled
- an output directory
- a Kraken2 16S database

The driver configuration documents the remaining thresholds, file names, and optional stages. The main pipeline entry point is:

```bash
bash src/argmapper/driver.sh
```

### Main outputs

Depending on enabled stages, ARGMapper writes barcode-level TSV summaries, ASV and sub-ARG summaries, filtered count tables, and an AMR/ASV–ARG heat map. The driver prints the output paths when processing finishes.

## Optional: PrimerPicker

PrimerPicker designs multiplex PCR primer pools from a FASTA file, scores primer dimerization, and ranks pools with simulated annealing. It is not needed if primer sequences are already available.

### Requirements

PrimerPicker requires Python 3.10 or newer and the `primer3_core` and `ntthal` executables:

```bash
conda create -n primer-picker -c bioconda -c conda-forge primer3 python=3.10
conda activate primer-picker
```

No Biopython, NumPy, joblib, or matplotlib installation is required by the CLI.

### Input and usage

Provide one target per FASTA record:

```text
>gene1
ATGCGTACGTAGCTAGCTAG...
>gene2
ATGAAACCCGGGTTTAAA...
```

From the repository root:

```bash
python PrimerPicker/primer_picker.py targets.fa \
  --outdir primer_picker_results \
  --seed 123
```

The default run generates up to 100 candidate pairs per target, scores dimers with `ntthal`, and performs 200 simulated-annealing runs. Run `python PrimerPicker/primer_picker.py --help` for all options.

The main output is `primer_picker_results/top-primer-sets.tsv`. The directory also contains candidate primers, FASTA sequences, dimer scores, and intermediate pickle files. The notebooks in `PrimerPicker/` preserve the original three-stage workflow; the CLI is the recommended entry point.

## Reproducibility and review

For reproducible analyses, record input files, software versions, command-line options, and random seeds. Primer candidates should be reviewed and experimentally validated before synthesis.

## Citation

Lan F, Saba J, Ross TD, Zhou Z, Krauska K, Anantharaman K, Landick R, Venturelli OS. Massively parallel single-cell sequencing of diverse microbial populations. *Nature Methods* 21, 228–235 (2024). https://doi.org/10.1038/s41592-023-02157-7

## License

This project is released under the MIT License. See [LICENSE](LICENSE).
