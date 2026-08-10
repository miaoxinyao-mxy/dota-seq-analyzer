# ARGMapper: AMR mapping for DoTA-Seq

ARGMapper was developed for analysis of single-cell targeted sequencing data based on the DoTA-Seq framework. It identifies bacterial taxa, detects antibiotic-resistance-gene (ARG) signals, and maps ARGs to the taxa and individual cells in a sample.

PrimerPicker is included as an optional utility for designing multiplex primer pools when validated primer sequences are not already available. It is not required for running ARGMapper.

## Repository layout

```text
src/argmapper/       ARGMapper pipeline modules and driver configuration
PrimerPicker/        Optional primer-design utility and original notebooks
pyproject.toml        Python package metadata
LICENSE               MIT license
```

## ARGMapper

### What ARGMapper does

ARGMapper processes paired-end sequencing reads through the following main stages:

1. Extracts reads associated with the targeted regions and sample barcodes.
2. Identifies the bacterial taxa represented in the sample.
3. Matches ARG-related reads to individual barcodes.
4. Summarizes taxonomic and ARG information for each cell.
5. Filters low-confidence or low-quality assignments.
6. Produces tables and visualizations showing the distribution of ARGs across bacterial taxa.

The main result is a taxa–ARG map that helps answer questions such as which bacterial taxa carry particular antibiotic-resistance genes and how those associations change between samples or time points.

### Installation

Create a Python environment with Python 3.13 or newer and install the package from the repository root:

```bash
python -m pip install .
```

Installing the Python package alone is not sufficient to run the complete pipeline. The workflow also requires external tools and databases, including Kraken2 and, when sub-ARG analysis is enabled, BLAST.

Before running the pipeline, edit `src/argmapper/driver.sh` to set:

- forward and reverse FASTQ paths
- the DoTA-Seq primer CSV
- the output directory
- taxonomic and sub-ARG database paths
- optional analysis settings

### Required inputs

The standard workflow expects:

- forward and reverse FASTQ files
- a DoTA-Seq primer CSV
- a taxonomic classification database
- a sub-ARG database FASTA when sub-ARG analysis is enabled

Run the configured pipeline with:

```bash
bash src/argmapper/driver.sh
```

### Outputs

Depending on the enabled analysis steps, ARGMapper writes:

- barcode-level summary tables
- bacterial taxonomic classifications
- ARG and sub-ARG summary tables
- filtered confidence and count tables
- taxa–ARG heat maps and related figures

The driver reports the output paths when processing finishes.

## Optional: PrimerPicker

PrimerPicker is an optional utility for designing multiplex primer pools when validated primer sequences are not already available.

### Installation

PrimerPicker requires Python 3.10 or newer and the `primer3_core` and `ntthal` executables:

```bash
conda create -n primer-picker -c bioconda -c conda-forge primer3 python=3.10
conda activate primer-picker
```

### Usage

Provide one target sequence per FASTA record and run:

```bash
python PrimerPicker/primer_picker.py targets.fa \
  --outdir primer_picker_results \
  --seed 123
```

Run `python PrimerPicker/primer_picker.py --help` for additional options.

### Outputs

The main output is `primer_picker_results/top-primer-sets.tsv`. The output directory also contains candidate primer sequences and dimer-scoring results. The notebooks in `PrimerPicker/` preserve the original workflow.

## Reproducibility and review

For reproducible analyses, record input files, software versions, command-line options, and random seeds. Primer candidates should be reviewed and experimentally validated before synthesis.

## Citation

ARGMapper was developed for analysis of single-cell targeted sequencing data based on the DoTA-Seq framework. For the underlying DoTA-Seq method, please cite:

Lan F, Saba J, Ross TD, Zhou Z, Krauska K, Anantharaman K, Landick R, Venturelli OS. Massively parallel single-cell sequencing of diverse microbial populations. *Nature Methods* 21, 228–235 (2024). https://doi.org/10.1038/s41592-023-02157-7

## License

This project is released under the MIT License. See [LICENSE](LICENSE).
