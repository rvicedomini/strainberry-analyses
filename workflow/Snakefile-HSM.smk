import os, glob, snakemake

# The main entry point of your workflow.
# After configuring, running snakemake -n in a clone of this repository should successfully execute a dry-run of the workflow.
#report: "report/workflow.rst"

# Allow users to fix the underlying OS via singularity.
#singularity: "docker://continuumio/miniconda3" # Allow users to fix the underlying OS via singularity.

# Rules requiring internet connection
localrules: hsm_dl_reads, hsm_dl_assembly, hsm_dl_references

snakemake.utils.validate(config, schema="schemas/config.schema.yaml")
sample = config["sample"]


def hsm_plots(sample):
    plot_files = []
    if nstrains == 2:
        plot_files = [
            f'results/{sample}/assembly_eval/refcoverage_Sthermophilus_NWC_2_1-flye.svg',
            f'results/{sample}/assembly_eval/refcoverage_Ldelbrueckii_NWC_2_2-flye.svg',
            f'results/{sample}/assembly_eval/refcoverage_Lhelveticus_NWC_2_3-flye.svg',
            f'results/{sample}/assembly_eval/refcoverage_Lhelveticus_NWC_2_4-flye.svg',
            f'results/{sample}/assembly_eval/refcoverage_Sthermophilus_NWC_2_1-canu.svg',
            f'results/{sample}/assembly_eval/refcoverage_Ldelbrueckii_NWC_2_2-canu.svg',
            f'results/{sample}/assembly_eval/refcoverage_Lhelveticus_NWC_2_3-canu.svg',
            f'results/{sample}/assembly_eval/refcoverage_Lhelveticus_NWC_2_4-canu.svg',
        ]
    return plot_files


rule all:
    input:
        # NWC2 - ONT resources
        expand('{basepath}/{filename}', basepath=config['read_basepath'], filename=config['read_filenames']),
        config['ref_csv'],
        # strainberry input: strain-oblivious assemblies & alignments
        f'results/{sample}/assemblies/flye.fa',
        f'results/{sample}/alignments/flye.bam',
        f'results/{sample}/assemblies/lathe-p1.fa',
        f'results/{sample}/alignments/lathe-p1.bam',
        # strainberry separation of lathe reference assembly
        f'results/{sample}/assemblies/sberry_lathe-p1_n{nstrains}_ctg.fa',
        f'results/{sample}/assemblies/sberry_lathe-p1_n{nstrains}_scf.fa',
        # assembly evaluation stats
        #f'results/{sample}/assembly_eval/flye.report.tsv',
        # reference coverage plots
        #nwc2_plots(sample),


# Common utilities (e.g., samtools faidx/index)
include: "rules/utils.smk"

# Retrieve resources and define read directory
include: "rules/hsm_resources.smk"

# Generation of input for Strainberry
include: "rules/hsm_assembly.smk"
include: "rules/alignment.smk"

# Strainberry rules
include: "rules/strainberry.smk"

# Evaluation rules against available references
include: "rules/hsm_evaluation.smk"

# Specific plots for Mock3 and Mock9
#include: "rules/nwc2_plots.smk"

