import os, glob, snakemake

# Rules requiring internet connection
localrules: hsm_dl_reads, hsm_dl_assembly, hsm_dl_references

sample = config['sample']

rule all:
    input:
        # NWC2 - ONT resources
        expand('{basepath}/{filename}', basepath=config['read_basepath'], filename=config['read_filenames']),
        config['ref_csv'],
        # Strainberry input: strain-oblivious assemblies & alignments
        f'results/{sample}/assemblies/flye.fa',
        f'results/{sample}/alignments/flye.bam',
        f'results/{sample}/assemblies/lathe-p1.fa',
        f'results/{sample}/alignments/lathe-p1.bam',
        # Strainberry separation of lathe reference assembly
        f'results/{sample}/assemblies/sberry_lathe-p1.fa',
        # Kraken2 classification of contigs
        f'results/{sample}/kraken2/lathe-p1.kraken2',
        # Lathe binning
        f'results/{sample}/binning/lathe-p1.depth.txt',
        f'results/{sample}/binning/lathe-p1.bins.tsv',
        f'results/{sample}/evaluation/lathe-p1.bin_stats.tsv',
        # Strainberry separation of lathe's best bins
        f'results/{sample}/evaluation/strainberry.bin_stats.tsv',
        f'results/{sample}/evaluation/strainberry_racon.bin_stats.tsv',
        f'results/{sample}/evaluation/strainberry_medaka.bin_stats.tsv',
        # V. atypica and E. eligens plots
        f'results/{sample}/evaluation/strainberry_lathe-p1.eeligens.depth.pdf',
        f'results/{sample}/evaluation/strainberry_lathe-p1.vatypica.depth.pdf',
        f'results/{sample}/evaluation/strainberry_medaka.vatypica.ACS-049-V-Sch6.ps',
        f'results/{sample}/evaluation/strainberry_medaka.vatypica.ACS-134-V-Col7a.ps',
        # pre/post separation plot
        f'results/{sample}/evaluation/lathe-p1_strainberry.barplot.pdf',
        f'results/{sample}/evaluation/lathe-p1_strainberry.barplot.svg',


# Common utilities (e.g., samtools faidx/index)
include: "rules/utils.smk"

# Retrieve resources and define read directory
include: "rules/hsm_resources.smk"

# Generation of input for Strainberry
include: "rules/hsm_assembly.smk"
include: "rules/alignment.smk"

# Canonical Strainberry separation
include: "rules/strainberry.smk"

# Bin-specific Strainberry separation and polishing
include: "rules/hsm_strainberry.smk"

# Classification, binning, evaluation, plots
include: "rules/hsm_evaluation.smk"


