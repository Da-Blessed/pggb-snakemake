import csv
from pathlib import Path


ROOT = Path(workflow.basedir).parent.resolve()
CHROMS = config["chromosomes"]
PARTITION_DIR = config["paths"]["partition_dir"]
REFERENCE_DIR = config["paths"]["reference_dir"]
MERGED_DIR = config["paths"]["merged_dir"]
OUTPUT_DIR = config["paths"]["output_dir"]
WORK_DIR = config["paths"]["work_dir"]
LOG_DIR = config["paths"]["log_dir"]
COMBINED_REF = f"{REFERENCE_DIR}/combined.fa"


def resolve(path):
    path = Path(path)
    return str(path if path.is_absolute() else ROOT / path)


with open(config["paths"]["samples"]) as handle:
    rows = list(csv.DictReader(handle, delimiter="\t"))

SAMPLES = {row["sample"]: resolve(row["fasta"]) for row in rows}
REFERENCES = [row["sample"] for row in rows if row["role"] == "reference"]
QUERIES = [row["sample"] for row in rows if row["role"] == "query"]

if len(REFERENCES) < 1 or len(QUERIES) < 1:
    raise ValueError("samples.tsv requires at least one reference and one query")


rule combined_reference:
    input:
        [SAMPLES[sample] for sample in REFERENCES]
    output:
        fasta=COMBINED_REF,
        fai=f"{COMBINED_REF}.fai"
    log:
        f"{LOG_DIR}/combined-reference.log"
    conda:
        "../envs/pggb.yaml"
    script:
        "../scripts/combine_fastas.py"


rule partition_sample:
    input:
        fasta=lambda wildcards: SAMPLES[wildcards.sample],
        reference=COMBINED_REF,
        reference_fai=f"{COMBINED_REF}.fai"
    output:
        chrom_fastas=[f"{PARTITION_DIR}/{chrom}/{{sample}}.fa" for chrom in CHROMS],
        assignment=f"{PARTITION_DIR}/assignments/{{sample}}.tsv"
    params:
        chromosomes=CHROMS,
        work_dir=lambda wildcards: f"{WORK_DIR}/partition/{wildcards.sample}",
        full_identity=config["partition"]["full_identity"],
        split_identity=config["partition"]["split_identity"],
        split_length=config["partition"]["split_length"],
        min_query_coverage=config["partition"]["min_query_coverage"]
    threads:
        config["resources"]["partition_threads"]
    resources:
        mem_mb=config["resources"]["partition_memory_mb"]
    log:
        f"{LOG_DIR}/partition/{{sample}}.log"
    wildcard_constraints:
        sample="|".join(QUERIES)
    conda:
        "../envs/pggb.yaml"
    script:
        "../scripts/partition_sample.py"


rule extract_reference:
    input:
        fasta=COMBINED_REF,
        fai=f"{COMBINED_REF}.fai"
    output:
        f"{PARTITION_DIR}/{{chrom}}/{{reference}}.fa"
    params:
        region=lambda wildcards: config["references"]["regions"][wildcards.reference].format(
            chrom=wildcards.chrom
        )
    wildcard_constraints:
        chrom="|".join(CHROMS),
        reference="|".join(REFERENCES)
    log:
        f"{LOG_DIR}/references/{{chrom}}.{{reference}}.log"
    conda:
        "../envs/pggb.yaml"
    shell:
        "mkdir -p $(dirname {output:q}) $(dirname {log:q}) && "
        "samtools faidx {input.fasta:q} {params.region:q} > {output:q} 2> {log:q}"


rule merge_chromosome:
    input:
        queries=lambda wildcards: [
            f"{PARTITION_DIR}/{wildcards.chrom}/{sample}.fa" for sample in QUERIES
        ],
        references=lambda wildcards: [
            f"{PARTITION_DIR}/{wildcards.chrom}/{sample}.fa" for sample in REFERENCES
        ]
    output:
        fasta=f"{MERGED_DIR}/{{chrom}}.merged.fa",
        fai=f"{MERGED_DIR}/{{chrom}}.merged.fa.fai",
        metadata=f"{MERGED_DIR}/{{chrom}}.metadata.json"
    log:
        f"{LOG_DIR}/merge/{{chrom}}.log"
    conda:
        "../envs/pggb.yaml"
    script:
        "../scripts/merge_chromosome.py"


rule pggb:
    input:
        fasta=f"{MERGED_DIR}/{{chrom}}.merged.fa",
        fai=f"{MERGED_DIR}/{{chrom}}.merged.fa.fai",
        metadata=f"{MERGED_DIR}/{{chrom}}.metadata.json"
    output:
        done=touch(f"{OUTPUT_DIR}/{{chrom}}/.pggb.done")
    params:
        output_dir=lambda wildcards: f"{OUTPUT_DIR}/{wildcards.chrom}",
        temp_dir=lambda wildcards: f"{WORK_DIR}/pggb/{wildcards.chrom}",
        identity=config["pggb"]["identity"],
        segment_length=config["pggb"]["segment_length"],
        references=config["pggb"]["vcf_references"],
        extra_args=config["pggb"].get("extra_args", [])
    threads:
        config["resources"]["pggb_threads"]
    resources:
        mem_mb=config["resources"]["pggb_memory_mb"]
    log:
        f"{LOG_DIR}/pggb/{{chrom}}.log"
    conda:
        "../envs/pggb.yaml"
    script:
        "../scripts/run_pggb.py"
