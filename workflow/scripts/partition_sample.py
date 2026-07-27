import gzip
import re
import shutil
import subprocess
from collections import defaultdict
from pathlib import Path


CHROM_PATTERN = re.compile(r"^(chr(?:[1-9]|1[0-9]|2[0-2]|X|Y|M|MT))(?:[_.].*)?$", re.I)


def normalize_chrom(name):
    match = CHROM_PATTERN.match(name.rsplit("#", 1)[-1])
    if not match:
        return None
    value = match.group(1).lower()
    if value in {"chrm", "chrmt"}:
        return "chrM"
    if value == "chrx":
        return "chrX"
    if value == "chry":
        return "chrY"
    return "chr" + value[3:]


def run(command, log, stdout=None):
    log.write("[cmd] " + " ".join(map(str, command)) + "\n")
    log.flush()
    subprocess.run(command, check=True, stdout=stdout or log, stderr=log)


work = Path(snakemake.params.work_dir)
work.mkdir(parents=True, exist_ok=True)
Path(snakemake.log[0]).parent.mkdir(parents=True, exist_ok=True)
source = Path(snakemake.input.fasta)
query = work / f"{snakemake.wildcards.sample}.fa"

with open(snakemake.log[0], "w") as log:
    if source.suffix == ".gz":
        with gzip.open(source, "rb") as src, query.open("wb") as dst:
            shutil.copyfileobj(src, dst)
    else:
        shutil.copyfile(source, query)

    run(["samtools", "faidx", str(query)], log)
    lengths = {}
    with open(str(query) + ".fai") as handle:
        for line in handle:
            name, length = line.split("\t")[:2]
            lengths[name] = int(length)

    full_paf = work / "full.paf"
    with full_paf.open("w") as paf:
        run([
            "wfmash", "-p", str(snakemake.params.full_identity), "-N",
            "-t", str(snakemake.threads), str(snakemake.input.reference), str(query),
        ], log, paf)

    assignment = {}
    scores = {}
    with full_paf.open() as handle:
        for line in handle:
            fields = line.rstrip().split("\t")
            if len(fields) < 12:
                continue
            qname = fields[0]
            coverage = (int(fields[3]) - int(fields[2])) / lengths[qname]
            chrom = normalize_chrom(fields[5])
            score = int(fields[10])
            if coverage >= float(snakemake.params.min_query_coverage) and chrom and score > scores.get(qname, -1):
                assignment[qname] = chrom
                scores[qname] = score

    unassigned = [name for name in lengths if name not in assignment]
    if unassigned:
        regions = work / "unassigned.txt"
        regions.write_text("\n".join(unassigned) + "\n")
        unassigned_fasta = work / "unassigned.fa"
        with unassigned_fasta.open("w") as fasta:
            run(["samtools", "faidx", str(query), "-r", str(regions)], log, fasta)
        split_paf = work / "split.paf"
        with split_paf.open("w") as paf:
            run([
                "wfmash", "-p", str(snakemake.params.split_identity),
                "-s", str(snakemake.params.split_length), "-t", str(snakemake.threads),
                str(snakemake.input.reference), str(unassigned_fasta),
            ], log, paf)
        aggregated = defaultdict(lambda: defaultdict(int))
        with split_paf.open() as handle:
            for line in handle:
                fields = line.rstrip().split("\t")
                if len(fields) < 11:
                    continue
                chrom = normalize_chrom(fields[5])
                if chrom:
                    aggregated[fields[0]][chrom] += int(fields[9])
        for name, values in aggregated.items():
            assignment[name] = max(values.items(), key=lambda item: item[1])[0]

    for name in lengths:
        assignment.setdefault(name, "unplaced")

    assignment_path = Path(snakemake.output.assignment)
    assignment_path.parent.mkdir(parents=True, exist_ok=True)
    with assignment_path.open("w") as handle:
        handle.write("sample\tcontig\tchromosome\n")
        for name, chrom in assignment.items():
            handle.write(f"{snakemake.wildcards.sample}\t{name}\t{chrom}\n")

    by_chrom = defaultdict(list)
    for name, chrom in assignment.items():
        by_chrom[chrom].append(name)
    for chrom, output_string in zip(snakemake.params.chromosomes, snakemake.output.chrom_fastas):
        output = Path(output_string)
        output.parent.mkdir(parents=True, exist_ok=True)
        names = by_chrom.get(chrom, [])
        if not names:
            output.write_text("")
            continue
        regions = work / f"{chrom}.txt"
        regions.write_text("\n".join(names) + "\n")
        with output.open("w") as fasta:
            run(["samtools", "faidx", str(query), "-r", str(regions)], log, fasta)
