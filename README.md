# DoTA-AMR: Mapping antibiotic resistance in single-cell DoTA-Seq data

DoTA-AMR links antimicrobial resistance genes to their bacterial hosts in single-cell DoTA-Seq data.

## Installation

Install DoTA-AMR and its reference databases:

```bash
git clone https://github.com/miaoxinyao-mxy/dota-amr.git
cd dota-amr
git lfs install
git lfs pull
conda env create -f environment.yml
tar -xzf database/dota-amr-taxonomy-db.tar.gz -C database
tar -xzf database/dota-amr-arg-db.tar.gz -C database
conda activate dota-seq-amr
```

## Inputs

- paired R1 and R2 FASTQ files
- a DoTA-Seq primer CSV file
- optionally, a custom AMR reference FASTA for matching detected sequences to known alleles

The taxonomic and AMR reference databases are included with the repository.

## Output

The primary results are `dota_amr_results.jsonl` and `reports/cell_amr_matrix.tsv`. The JSONL contains one structured record per cell; the TSV contains one row per cell and one column per AMR target, together with taxonomic and quality-control information.

```bash
python src/dota_amr/export_results.py \
  --input_tsv reports/cell_amr_matrix.tsv \
  --primers_file primers.csv \
  --output_jsonl dota_amr_results.jsonl
```

Intermediate files are written to `tmp/`, report tables to `reports/`, and figures to `figures/`.

## Optional: PrimerPicker

Provide one target sequence per FASTA record and run:

```bash
python PrimerPicker/primer_picker.py targets.fa \
  --outdir primer_picker_results \
  --seed 123
```

The main output is `primer_picker_results/top-primer-sets.tsv`. Run `python PrimerPicker/primer_picker.py --help` for additional options.

## Citation

DoTA-AMR was developed for analysis of single-cell targeted sequencing data based on the DoTA-Seq framework. For the underlying DoTA-Seq method, please cite:

Lan F, Saba J, Ross TD, Zhou Z, Krauska K, Anantharaman K, Landick R, Venturelli OS. Massively parallel single-cell sequencing of diverse microbial populations. *Nature Methods* 21, 228–235 (2024). https://doi.org/10.1038/s41592-023-02157-7

## License

This project is released under the MIT License. See [LICENSE](LICENSE).
