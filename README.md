# DoTA-Seq tools: ARGMapper and PrimerPicker

This repository contains two companion tools for the DoTA-Seq workflow:

- **ARGMapper** processes paired-end DoTA-Seq reads and links antibiotic-resistance-gene (ARG) signal to bacterial taxonomic classifications.
- **PrimerPicker** designs multiplex PCR primer pools from a FASTA file, scores primer dimerization, and ranks pools with simulated annealing.

The tools are useful together: PrimerPicker can produce candidate pools for a DoTA-Seq experiment, and ARGMapper can analyze the resulting sequencing data at one or more time points to investigate ARG carriage and horizontal gene transfer.

## Repository layout

```text
src/argmapper/       ARGMapper pipeline modules and driver configuration
PrimerPicker/        PrimerPicker CLI and original workflow notebooks
pyproject.toml        Python package metadata
LICENSE               MIT license
```

## PrimerPicker

### Requirements

PrimerPicker requires Python 3.10 or newer and the `primer3_core` and `ntthal` executables. The simplest installation is through Conda:

```bash
conda create -n primer-picker -c bioconda -c conda-forge primer3 python=3.10
conda activate primer-picker
```

No Biopython, NumPy, joblib, or matplotlib installation is required by the CLI.

### Input

Provide one target per FASTA record. The record name becomes the target name:

```text
>gene1
ATGCGTACGTAGCTAGCTAG...
>gene2
ATGAAACCCGGGTTTAAA...
```

### Run

From the repository root:

```bash
python PrimerPicker/primer_picker.py targets.fa \
  --outdir primer_picker_results \
  --seed 123
```

The default run generates up to 100 candidate pairs per target, scores dimers with `ntthal`, and performs 200 simulated-annealing runs. For a quick test:

```bash
python PrimerPicker/primer_picker.py targets.fa \
  --num-primers 10 \
  --anneal-iterations 5 \
  --anneal-minutes 0.01 \
  --dimer-jobs 1 \
  --anneal-jobs 1
```

Useful options include `--product-size-range 400-500`, repeated `--exclude TEXT`, `--no-adapters`, `--no-16s`, and explicit executable paths with `--primer3-core PATH` and `--ntthal PATH`. Run `python PrimerPicker/primer_picker.py --help` for the complete list.

### Results

The output directory contains:

- `top-primer-sets.tsv`: ranked primer pools for experimental review
- `primers.tsv` and `primers.fasta`: all generated candidates
- `dimerization-deltaGs.tsv`: calculated dimer scores
- pickle files used to resume or inspect intermediate results

The notebooks in `PrimerPicker/` preserve the original Primer3, dimer-scoring, and simulated-annealing workflow. They are optional; the CLI is the recommended entry point.

## ARGMapper

ARGMapper is intended for paired-end, short-read bacterial samples generated as part of the DoTA-Seq workflow. It uses 16S taxonomic classification, barcode-level processing, ASV information, ARG filtering, optional sub-ARG analysis, and statistical classification to produce per-cell summaries and ASV–ARG visualizations.

### Installation

Create a Python environment with Python 3.13 or newer and install the package from the repository root:

```bash
python -m pip install .
```

The full pipeline also depends on external DoTA-Seq tools and databases, including Kraken2 and (when sub-ARG analysis is enabled) BLAST. Configure their paths and analysis parameters in `src/argmapper/driver.sh` before running the workflow.

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

Depending on enabled stages, ARGMapper writes barcode-level TSV summaries, ASV and sub-ARG summaries, filtered count tables, and an ASV–ARG heat map. The driver prints the output paths when processing finishes.

## Reproducibility and review

Primer selection is computational support for experimental design; inspect the ranked pools and validate them experimentally before synthesis. For reproducible PrimerPicker runs, record the input FASTA, executable versions, command-line options, and random seed.

## Citation

Lan F, Saba J, Ross TD, Zhou Z, Krauska K, Anantharaman K, Landick R, Venturelli OS. Massively parallel single-cell sequencing of diverse microbial populations. *Nature Methods* 21, 228–235 (2024). https://doi.org/10.1038/s41592-023-02157-7

## License

This project is released under the MIT License. See [LICENSE](LICENSE).
