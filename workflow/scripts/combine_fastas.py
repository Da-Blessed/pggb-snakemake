import gzip
import shutil
import subprocess
from pathlib import Path


output = Path(snakemake.output.fasta)
output.parent.mkdir(parents=True, exist_ok=True)
Path(snakemake.log[0]).parent.mkdir(parents=True, exist_ok=True)
with output.open("wb") as destination:
    for path_string in snakemake.input:
        path = Path(path_string)
        opener = gzip.open if path.suffix == ".gz" else open
        with opener(path, "rb") as source:
            shutil.copyfileobj(source, destination)

with open(snakemake.log[0], "w") as log:
    subprocess.run(
        ["samtools", "faidx", str(output)],
        check=True,
        stdout=log,
        stderr=subprocess.STDOUT,
    )
