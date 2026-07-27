import json
import shlex
import subprocess
from pathlib import Path


metadata = json.loads(Path(snakemake.input.metadata).read_text())
output_dir = Path(snakemake.params.output_dir)
temp_dir = Path(snakemake.params.temp_dir)
output_dir.mkdir(parents=True, exist_ok=True)
temp_dir.mkdir(parents=True, exist_ok=True)
Path(snakemake.log[0]).parent.mkdir(parents=True, exist_ok=True)

command = [
    "pggb",
    "-i", str(snakemake.input.fasta),
    "-o", str(output_dir),
    "-n", str(metadata["sample_count"]),
    "-t", str(snakemake.threads),
    "-p", str(snakemake.params.identity),
    "-s", str(snakemake.params.segment_length),
    "-D", str(temp_dir),
]
references = list(snakemake.params.references)
if references:
    command.extend(["-V", ",".join(references)])
command.extend(map(str, snakemake.params.extra_args))

with open(snakemake.log[0], "w") as log:
    log.write("[cmd] " + shlex.join(command) + "\n")
    log.flush()
    subprocess.run(command, check=True, stdout=log, stderr=subprocess.STDOUT)
