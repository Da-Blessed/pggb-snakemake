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

### FASTA sequence naming

All input FASTA sequence names must already follow the PanSN convention:

```text
sample#haplotype#contig
```

For example:

```text
>HG00438#1#chr1
>HG00438#2#chr1
>CHM13#1#chr1
```

Sequence names must be unique across all input FASTAs. This workflow preserves the existing names and does not automatically convert non-PanSN headers. Use only the first whitespace-delimited field of a FASTA header as the sequence name, and ensure that its `sample` and `haplotype` fields correctly identify the path before running the workflow.

With the default configuration, each partition or PGGB job uses 8 threads. Running with 32 cores allows up to four jobs concurrently.

## Run

```bash
conda env create -f environment.yaml
conda activate pggb-snakemake

snakemake \
  --cores 32 \
  --use-conda \
  --conda-prefix .conda/rules \
  --printshellcmds
```

Inspect the DAG without running jobs:

```bash
snakemake \
  --dry-run \
  --cores 32 \
  --use-conda \
  --conda-prefix .conda/rules
```

Run these commands from the repository root. Adjust `--cores` to the resources available on the machine. The `--conda-prefix` option keeps rule environments under `.conda/rules` in the repository.

To keep Snakemake metadata and locks in a separate run directory, use an absolute
Snakefile path together with `--directory`. The workflow still loads
`config/config.yaml` from the cloned repository:

```bash
REPO=/path/to/pggb-snakemake
RUN_DIR=/path/to/pggb-run

mkdir -p "$RUN_DIR"
snakemake \
  --snakefile "$REPO/Snakefile" \
  --directory "$RUN_DIR" \
  --cores 32 \
  --use-conda \
  --conda-prefix "$REPO/.conda/rules" \
  --printshellcmds
```

If the user home directory is read-only, set a writable Snakemake cache location before running:

```bash
export XDG_CACHE_HOME="$PWD/.cache"
```

Per-sample partition logs and per-chromosome PGGB logs are written below `logs/`. Intermediate partition and merged FASTA files are written below `work/`.
