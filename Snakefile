configfile: workflow.basedir + "/config/config.yaml"

include: "workflow/rules/pggb.smk"


rule all:
    input:
        expand(f"{OUTPUT_DIR}/{{chrom}}/.pggb.done", chrom=CHROMS)
    default_target:
        True
