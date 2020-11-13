from snakemake.utils import validate
#import pandas as pd

# this container defines the underlying OS for each job when using the workflow
# with --use-conda --use-singularity
singularity: "docker://continuumio/miniconda3"

##### load config and sample sheets #####
configfile: "config/config.yaml"

validate(config, schema="../schemas/config.schema.yaml")
sample = config["sample"]
nstrains = int(config["nstrains"])

rule samtools_faidx:
    input:  '{something}.f{asta}'
    output: '{something}.f{asta}.fai'
    conda:  '../envs/common.yaml'
    shell:  'samtools faidx {input}'


rule samtools_index:
    input:  '{something}.bam'
    output: '{something}.bam.bai'
    conda:  '../envs/common.yaml'
    shell:  'samtools index -@ 4 {input}'




