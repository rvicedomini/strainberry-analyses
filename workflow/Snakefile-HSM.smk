import os, glob, snakemake

# The main entry point of your workflow.
# After configuring, running snakemake -n in a clone of this repository should successfully execute a dry-run of the workflow.
#report: "report/workflow.rst"

# Allow users to fix the underlying OS via singularity.
#singularity: "docker://continuumio/miniconda3" # Allow users to fix the underlying OS via singularity.

# Rules requiring internet connection
localrules: hsm_dl_reads, hsm_dl_assembly, hsm_dl_references

snakemake.utils.validate(config, schema="schemas/config.schema.yaml")
sample = config['sample']
nstrains = int(config['nstrains'])

#def hsm_plots(sample):
#    plot_files = []
#    if nstrains == 2:
#        plot_files = [
#            f'results/{sample}/assembly_eval/refcoverage_Sthermophilus_NWC_2_1-flye.svg',
#            f'results/{sample}/assembly_eval/refcoverage_Ldelbrueckii_NWC_2_2-flye.svg',
#            f'results/{sample}/assembly_eval/refcoverage_Lhelveticus_NWC_2_3-flye.svg',
#            f'results/{sample}/assembly_eval/refcoverage_Lhelveticus_NWC_2_4-flye.svg',
#            f'results/{sample}/assembly_eval/refcoverage_Sthermophilus_NWC_2_1-canu.svg',
#            f'results/{sample}/assembly_eval/refcoverage_Ldelbrueckii_NWC_2_2-canu.svg',
#            f'results/{sample}/assembly_eval/refcoverage_Lhelveticus_NWC_2_3-canu.svg',
#            f'results/{sample}/assembly_eval/refcoverage_Lhelveticus_NWC_2_4-canu.svg',
#        ]
#    return plot_files


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
        # strainberry polished contig assembly
        f'results/{sample}/assemblies/medaka/sberry_lathe-p1_n{nstrains}_ctg/contigs',
        f'results/{sample}/assemblies/sberry_lathe-p1_n{nstrains}_ctg.medaka.fa',
        # kraken2 classification of contigs
        f'results/{sample}/kraken2/lathe-p1.kraken2',
#        f'results/{sample}/kraken2/sberry_lathe-p1_n{nstrains}_ctg.medaka.kraken2',
        # binning
        f'results/{sample}/binning/lathe-p1.depth.txt',
#        f'results/{sample}/binning/sberry_lathe-p1_n{nstrains}_ctg.medaka.depth.txt',
        # TODO best bins
        #f'results/{sample}/assembly_eval/lathe-p1.best.tsv',
        #f'results/{sample}/assembly_eval/sberry_lathe-p1_n{nstrains}_ctg.best.tsv',
        #f'results/{sample}/assembly_eval/sberry_lathe-p1_n{nstrains}_ctg.polished.best.tsv',
        # TODO pre/post separation plot
        #f'results/{sample}/assembly_eval/sberry_lathe-p1_n{nstrains}_ctg.barplot.pdf',



# Common utilities (e.g., samtools faidx/index)
include: "rules/utils.smk"

# Retrieve resources and define read directory
include: "rules/hsm_resources.smk"

# Generation of input for Strainberry
include: "rules/hsm_assembly.smk"
include: "rules/alignment.smk"

# Strainberry rules
include: "rules/strainberry.smk"

# Polishing
include: 'rules/hsm_polishing.smk'

# Evaluation rules of generated assemblies (classification, binning, etc.)
include: "rules/hsm_evaluation.smk"

# Specific plots for HSM
#include: "rules/hsm_plots.smk"

