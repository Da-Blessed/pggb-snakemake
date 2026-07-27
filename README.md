# PGGB Snakemake workflow

Reproducible chromosome-partitioned workflow for PGGB 0.7.4. It assigns query contigs to CHM13/GRCh38 chromosomes with wfmash, adds the two reference paths, merges chromosome FASTAs, and runs PGGB independently for each chromosome.

## Requirements

- Linux x86-64
- Conda/Miniforge
- Sufficient local resources
- Input FASTA files are external data and are not committed to Git

PGGB and its core graph dependencies are pinned in `workflow/envs/pggb.yaml`.

## Configure

Edit `config/samples.tsv` and `config/config.yaml`. Relative FASTA paths are resolved from the parent directory of this workflow repository; absolute paths are also accepted.

`role` must be either `reference` or `query`. Reference sequence regions are configured under `references.regions`.

With the default configuration, each partition or PGGB job uses 8 threads. Running with 32 cores allows up to four jobs concurrently.

## Run

```bash
conda env create -f environment.yaml
conda activate pggb-snakemake

snakemake \
  --cores 32 \
  --use-conda \
  --conda-frontend conda \
  --printshellcmds
```

Inspect the DAG without running jobs:

```bash
snakemake --dry-run --cores 32 --use-conda
```

Per-sample partition logs and per-chromosome PGGB logs are written below `logs/`. Intermediate partition and merged FASTA files are written below `work/`.
