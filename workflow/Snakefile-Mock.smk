import os, glob, snakemake

# Rules requiring internet connection
localrules: mock_dl_references, mock_dl_data, mock_dl_barcodes

sample = config["sample"]

def mock_plots(sample):
    plot_files = []
    if sample == 'mock3':
        plot_files = [
            f"results/{sample}/evaluation/{sample}.ani.barplot.pdf",
            f"results/{sample}/evaluation/{sample}.dupratio.barplot.pdf",
            f"results/{sample}/evaluation/{sample}_circos.svg",
            ]
    elif sample == 'mock9':
        plot_files = [
            f"results/{sample}/evaluation/{sample}_n2.ani.barplot.pdf",
            f"results/{sample}/evaluation/{sample}_n2.dupratio.barplot.pdf",
            f"results/{sample}/evaluation/{sample}_n2_circos.svg",
            f"results/{sample}/evaluation/{sample}_n3.ani.barplot.pdf",
            f"results/{sample}/evaluation/{sample}_n3.dupratio.barplot.pdf",
            f"results/{sample}/evaluation/{sample}_n3_circos.svg",
            ]
    return plot_files


rule all:
    input:
        # Mock3 & Mock9 resources
        expand('{basepath}/{filename}', basepath=config['read_basepath'], filename=config['read_filenames']),
        config['ref_csv'],
        # strainberry input: strain-oblivious assemblies & alignments
        f'results/{sample}/assemblies/flye.fa',
        f'results/{sample}/assemblies/flye.fa.fai',
        f'results/{sample}/alignments/flye.bam',
        f'results/{sample}/assemblies/canu.fa',
        f'results/{sample}/assemblies/canu.fa.fai',
        f'results/{sample}/alignments/canu.bam',
        # strainberry assemblies
        f'results/{sample}/assemblies/sberry_flye.fa',
        f'results/{sample}/assemblies/sberry_canu.fa',
        # assembly evaluation stats
        f'results/{sample}/evaluation/flye.report.tsv',
        f'results/{sample}/evaluation/sberry_flye.report.tsv',
        f'results/{sample}/evaluation/canu.report.tsv',
        f'results/{sample}/evaluation/sberry_canu.report.tsv',
        # assembly evaluation barplots & circos
        mock_plots(sample),

# Common utilities (e.g., samtools faidx/index)
include: "rules/utils.smk"

# Retrieve resources and define read directory
include: "rules/mock_resources.smk"

# Generation of input for Strainberry
include: "rules/assembly.smk"
include: "rules/alignment.smk"

# Strainberry rules
include: "rules/strainberry.smk"

# Evaluation rules against available references
include: "rules/mock_evaluation.smk"

# Specific plots for Mock3 and Mock9
include: "rules/mock3_plots.smk"
include: "rules/mock9_n2_plots.smk"
include: "rules/mock9_n3_plots.smk"

