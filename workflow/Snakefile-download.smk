# The main entry point of your workflow.
# After configuring, running snakemake -n in a clone of this repository should successfully execute a dry-run of the workflow.
#report: "report/workflow.rst"

localrules: download_mock, download_nwc2_pacbio, download_nwc2_ont, download_hsm

# Common rules, such as fasta/bam indexing, etc.
#include: "rules/common.smk"

# Allow users to fix the underlying OS via singularity.
#singularity: "docker://continuumio/miniconda3" # Allow users to fix the underlying OS via singularity.


rule all:
    input:
        "resources/mock9/source/MM_10Plex_42MB_SEQUENCING_DATA_EXPRESS2.0.tar.gz",
        "resources/nwc2_pacbio/source/SRR7585901.fastq.gz",
        "resources/nwc2_ont/source/SRR7585900.fastq.gz",


rule download_mock:
    input: ancient("config/mock9.yaml")
    output: "resources/mock9/source/MM_10Plex_42MB_SEQUENCING_DATA_EXPRESS2.0.tar.gz"
    shell: 
        """
        curl --create-dirs -L -o resources/mock/source/MM_10Plex_42MB_SEQUENCING_DATA_EXPRESS2.0.tar.gz https://downloads.pacbcloud.com/public/dataset/microbial_multiplex_dataset_release_SMRT_Link_v6.0.0_with_Express_2.0/MM_10Plex_42MB_SEQUENCING_DATA_EXPRESS2.0.tar.gz
        """

rule download_nwc2_pacbio:
    input: ancient("config/nwc2_pacbio.yaml")
    output: "resources/nwc2_pacbio/source/SRR7585901.fastq.gz"
    conda: "envs/preprocess.yaml"
    shell:
        """
        fasterq-dump -O resources/nwc2_pacbio/source -o SRR7585901.fastq --min-read-len 5000 SRR7585901 && gzip resources/nwc2_pacbio/source/SRR7585901.fastq
        """

rule download_nwc2_ont:
    input: ancient("config/nwc2_ont.yaml")
    output: "resources/nwc2_ont/source/SRR7585900.fastq.gz"
    conda: "envs/preprocess.yaml"
    shell:
        """
        fasterq-dump -O resources/nwc2_ont/source -o SRR7585900.fastq --min-read-len 10000 SRR7585900 && gzip resources/nwc2_ont/source/SRR7585900.fastq
        """



