import json
import shutil
import subprocess
from pathlib import Path


inputs = [Path(path) for path in list(snakemake.input.queries) + list(snakemake.input.references)]
nonempty = [path for path in inputs if path.stat().st_size > 0]
if not nonempty:
    raise RuntimeError(f"No non-empty FASTA inputs for {snakemake.wildcards.chrom}")

output = Path(snakemake.output.fasta)
output.parent.mkdir(parents=True, exist_ok=True)
Path(snakemake.log[0]).parent.mkdir(parents=True, exist_ok=True)
with output.open("wb") as destination:
    for path in nonempty:
        with path.open("rb") as source:
            shutil.copyfileobj(source, destination)

with open(snakemake.log[0], "w") as log:
    subprocess.run(
        ["samtools", "faidx", str(output)],
        check=True,
        stdout=log,
        stderr=subprocess.STDOUT,
    )

Path(snakemake.output.metadata).write_text(json.dumps({
    "chromosome": str(snakemake.wildcards.chrom),
    "sample_count": len(nonempty),
    "samples": [path.stem for path in nonempty],
}, indent=2) + "\n")
