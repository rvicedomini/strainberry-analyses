# The main entry point of your workflow.
# After configuring, running snakemake -n in a clone of this repository should successfully execute a dry-run of the workflow.
#report: "report/workflow.rst"

# Common rules, such as fasta/bam indexing, etc.
include: "rules/common.smk"

# Allow users to fix the underlying OS via singularity.
#singularity: "docker://continuumio/miniconda3" # Allow users to fix the underlying OS via singularity.


rule all:
    input:
        # input data
        # strain-oblivious assemblies
        f'results/{sample}/assemblies/canu.fa',
        f'results/{sample}/assemblies/flye.fa',
        # strainberry assemblies
        f'results/{sample}/assemblies/sberry_canu_n{nstrains}.contigs.fa',
        f'results/{sample}/assemblies/sberry_canu_n{nstrains}.scaffolds.fa',
        f'results/{sample}/assemblies/sberry_flye_n{nstrains}.contigs.fa',
        f'results/{sample}/assemblies/sberry_flye_n{nstrains}.scaffolds.fa',
        # assembly evaluation
        f'results/{sample}/assembly_eval/canu.report.tsv',
        f'results/{sample}/assembly_eval/flye.report.tsv',
        f'results/{sample}/assembly_eval/sberry_canu_n{nstrains}.contigs.report.tsv',
        f'results/{sample}/assembly_eval/sberry_canu_n{nstrains}.scaffolds.report.tsv',
        f'results/{sample}/assembly_eval/sberry_flye_n{nstrains}.contigs.report.tsv',
        f'results/{sample}/assembly_eval/sberry_flye_n{nstrains}.scaffolds.report.tsv',


# Input generation from raw data
#include: "rules/assembly.smk"
#include: "rules/alignment.smk"

# Strainberry rules
include: "rules/sberry.smk"

# Evaluation rules against available references
include: "rules/assembly_eval.smk"
#include: "rules/plots.smk"

