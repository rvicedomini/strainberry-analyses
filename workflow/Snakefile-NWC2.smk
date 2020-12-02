import os, glob, snakemake

# The main entry point of your workflow.
# After configuring, running snakemake -n in a clone of this repository should successfully execute a dry-run of the workflow.
#report: "report/workflow.rst"

# Allow users to fix the underlying OS via singularity.
#singularity: "docker://continuumio/miniconda3" # Allow users to fix the underlying OS via singularity.

# Rules requiring internet connection
localrules: nwc2_dl_references, nwc2_ont_dl_reads, nwc2_pacbio_dl_reads

#snakemake.utils.validate(config, schema="schemas/config.schema.yaml")
sample = config["sample"]
nstrains = int(config["nstrains"])


def nwc2_plots(sample):
    plot_files = []
    if nstrains == 2:
        plot_files = [
            f'results/{sample}/evaluation/refcoverage_Sthermophilus_NWC_2_1-flye.svg',
            f'results/{sample}/evaluation/refcoverage_Ldelbrueckii_NWC_2_2-flye.svg',
            f'results/{sample}/evaluation/refcoverage_Lhelveticus_NWC_2_3-flye.svg',
            f'results/{sample}/evaluation/refcoverage_Lhelveticus_NWC_2_4-flye.svg',
            f'results/{sample}/evaluation/refcoverage_Sthermophilus_NWC_2_1-canu.svg',
            f'results/{sample}/evaluation/refcoverage_Ldelbrueckii_NWC_2_2-canu.svg',
            f'results/{sample}/evaluation/refcoverage_Lhelveticus_NWC_2_3-canu.svg',
            f'results/{sample}/evaluation/refcoverage_Lhelveticus_NWC_2_4-canu.svg',
        ]
    return plot_files


rule all:
    input:
        # NWC2 - reads and references
        expand('{basepath}/{filename}', basepath=config['read_basepath'], filename=config['read_filenames']),
        config['ref_csv'],
        # strainberry input: strain-oblivious assemblies & alignments
        f'results/{sample}/assemblies/flye.fa',
        f'results/{sample}/alignments/flye.bam',
        f'results/{sample}/assemblies/canu.fa',
        f'results/{sample}/alignments/canu.bam',
        # strainberry assemblies
        f'results/{sample}/assemblies/sberry_flye_n{nstrains}_ctg.fa',
        f'results/{sample}/assemblies/sberry_flye_n{nstrains}_scf.fa',
        f'results/{sample}/assemblies/sberry_canu_n{nstrains}_ctg.fa',
        f'results/{sample}/assemblies/sberry_canu_n{nstrains}_scf.fa',
        # assembly evaluation stats
        f'results/{sample}/evaluation/flye.report.tsv',
        f'results/{sample}/evaluation/sberry_flye_n{nstrains}_ctg.report.tsv',
        f'results/{sample}/evaluation/sberry_flye_n{nstrains}_scf.report.tsv',
        f'results/{sample}/evaluation/canu.report.tsv',
        f'results/{sample}/evaluation/sberry_canu_n{nstrains}_ctg.report.tsv',
        f'results/{sample}/evaluation/sberry_canu_n{nstrains}_scf.report.tsv',
        # checkm evaluation
        f'results/{sample}/evaluation/flye.checkm.tsv',
        f'results/{sample}/evaluation/sberry_flye_n{nstrains}_ctg.checkm.tsv',
        f'results/{sample}/evaluation/sberry_flye_n{nstrains}_scf.checkm.tsv',
        f'results/{sample}/evaluation/canu.checkm.tsv',
        f'results/{sample}/evaluation/sberry_canu_n{nstrains}_ctg.checkm.tsv',
        f'results/{sample}/evaluation/sberry_canu_n{nstrains}_scf.checkm.tsv',
        # reference coverage plots
        nwc2_plots(sample),


# Common utilities (e.g., samtools faidx/index)
include: "rules/utils.smk"

# Retrieve resources and define read directory
include: "rules/nwc2_resources.smk"

# Generation of input for Strainberry
include: "rules/nwc2_assembly.smk"
include: "rules/alignment.smk"

# Strainberry rules
include: "rules/strainberry.smk"

# Evaluation rules against available references
include: "rules/nwc2_evaluation.smk"

# Specific plots for NWC2
include: "rules/nwc2_plots.smk"

