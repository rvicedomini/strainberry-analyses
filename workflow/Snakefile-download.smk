# The main entry point of your workflow.
# After configuring, running snakemake -n in a clone of this repository should successfully execute a dry-run of the workflow.
#report: "report/workflow.rst"

localrules: mock9_data, nwc2_pacbio_reads, nwc2_pacbio_references, nwc2_ont_reads, nwc2_ont_references

# Common rules, such as fasta/bam indexing, etc.
#include: "rules/common.smk"

# Allow users to fix the underlying OS via singularity.
#singularity: "docker://continuumio/miniconda3" # Allow users to fix the underlying OS via singularity.


rule all:
    input:
        # Mock datasets
        #"resources/mock9/source/MM_10Plex_42MB_SEQUENCING_DATA_EXPRESS2.0.tar.gz",
        #"resources/mock9/source/Sequel_16_Barcodes_v3.fasta",
        "resources/mock9/reads/lima/m54081_181221_163846.bc1001_BAK8A_OA--bc1001_BAK8A_OA.fq.gz",
        "resources/mock9/reads/lima/m54081_181221_163846.bc1002_BAK8A_OA--bc1002_BAK8A_OA.fq.gz",
        "resources/mock9/reads/lima/m54081_181221_163846.bc1009_BAK8A_OA--bc1009_BAK8A_OA.fq.gz",
        "resources/mock9/reads/lima/m54081_181221_163846.bc1010_BAK8A_OA--bc1010_BAK8A_OA.fq.gz",
        "resources/mock9/reads/lima/m54081_181221_163846.bc1012_BAK8A_OA--bc1012_BAK8A_OA.fq.gz",
        "resources/mock9/reads/lima/m54081_181221_163846.bc1015_BAK8B_OA--bc1015_BAK8B_OA.fq.gz",
        "resources/mock9/reads/lima/m54081_181221_163846.bc1016_BAK8B_OA--bc1016_BAK8B_OA.fq.gz",
        "resources/mock9/reads/lima/m54081_181221_163846.bc1018_BAK8B_OA--bc1018_BAK8B_OA.fq.gz",
        "resources/mock9/reads/lima/m54081_181221_163846.bc1019_BAK8B_OA--bc1019_BAK8B_OA.fq.gz",
        "resources/mock9/reads/lima/m54081_181221_163846.bc1022_BAK8B_OA--bc1022_BAK8B_OA.fq.gz"
#        # NWC2 - PacBio
#        "resources/nwc2_pacbio/source/SRR7585901.fastq.gz",
#        # NWC2 - Nanopore
#        "resources/nwc2_ont/source/SRR7585900.fastq.gz",


#####################
### MOCK9 DATASET ###
#####################

rule mock9_data:
    input: ancient("config/mock9.yaml")
    output: "resources/mock9/source/MM_10Plex_42MB_SEQUENCING_DATA_EXPRESS2.0.tar.gz"
    shell: 
        """
        curl --create-dirs -L -o resources/mock9/source/MM_10Plex_42MB_SEQUENCING_DATA_EXPRESS2.0.tar.gz \
                https://downloads.pacbcloud.com/public/dataset/microbial_multiplex_dataset_release_SMRT_Link_v6.0.0_with_Express_2.0/MM_10Plex_42MB_SEQUENCING_DATA_EXPRESS2.0.tar.gz
        """

rule mock9_barcodes:
    input:  ancient("config/mock9.yaml")
    output: "resources/mock9/source/Sequel_16_Barcodes_v3.fasta"
    shell:  
        """
        curl --create-dirs -s -L -o resources/mock9/source/Sequel_16_Barcodes_v3.zip https://www.pacb.com/wp-content/uploads/Sequel_16_Barcodes_v3.zip \
                && unzip -d resources/mock9/source/ resources/mock9/source/Sequel_16_Barcodes_v3.zip \
                && rm resources/mock9/source/Sequel_16_Barcodes_v3.zip
        """

rule mock9_demultiplex:
    input: 
        ancient("resources/mock9/source/MM_10Plex_42MB_SEQUENCING_DATA_EXPRESS2.0.tar.gz"),
        ancient("resources/mock9/source/Sequel_16_Barcodes_v3.fasta")
    output: 
        "resources/mock9/reads/lima/m54081_181221_163846.bc1001_BAK8A_OA--bc1001_BAK8A_OA.bam",
        "resources/mock9/reads/lima/m54081_181221_163846.bc1002_BAK8A_OA--bc1002_BAK8A_OA.bam",
        "resources/mock9/reads/lima/m54081_181221_163846.bc1009_BAK8A_OA--bc1009_BAK8A_OA.bam",
        "resources/mock9/reads/lima/m54081_181221_163846.bc1010_BAK8A_OA--bc1010_BAK8A_OA.bam",
        "resources/mock9/reads/lima/m54081_181221_163846.bc1012_BAK8A_OA--bc1012_BAK8A_OA.bam",
        "resources/mock9/reads/lima/m54081_181221_163846.bc1015_BAK8B_OA--bc1015_BAK8B_OA.bam",
        "resources/mock9/reads/lima/m54081_181221_163846.bc1016_BAK8B_OA--bc1016_BAK8B_OA.bam",
        "resources/mock9/reads/lima/m54081_181221_163846.bc1018_BAK8B_OA--bc1018_BAK8B_OA.bam",
        "resources/mock9/reads/lima/m54081_181221_163846.bc1019_BAK8B_OA--bc1019_BAK8B_OA.bam",
        "resources/mock9/reads/lima/m54081_181221_163846.bc1022_BAK8B_OA--bc1022_BAK8B_OA.bam"
    conda: 'envs/pacbio.yaml'
    threads: 
        workflow.cores
    shell:
        """
        tar -xf resources/mock9/source/MM_10Plex_42MB_SEQUENCING_DATA_EXPRESS2.0.tar.gz -C resources/mock9/source \
                && mkdir -p resources/mock9/reads/lima \
                && lima --num-threads {threads} --same --split-named \
                    resources/mock9/source/MM_10Plex_42MB_SEQUENCING_DATA_EXPRESS2.0/m54081_181221_163846.subreads.bam \
                    resources/mock9/source/Sequel_16_Barcodes_v3.fasta \
                    resources/mock9/reads/lima/m54081_181221_163846.bam
        """

rule mock9_fastq:
    input:
        "resources/mock9/reads/lima/m54081_181221_163846.bc1001_BAK8A_OA--bc1001_BAK8A_OA.bam",
        "resources/mock9/reads/lima/m54081_181221_163846.bc1002_BAK8A_OA--bc1002_BAK8A_OA.bam",
        "resources/mock9/reads/lima/m54081_181221_163846.bc1009_BAK8A_OA--bc1009_BAK8A_OA.bam",
        "resources/mock9/reads/lima/m54081_181221_163846.bc1010_BAK8A_OA--bc1010_BAK8A_OA.bam",
        "resources/mock9/reads/lima/m54081_181221_163846.bc1012_BAK8A_OA--bc1012_BAK8A_OA.bam",
        "resources/mock9/reads/lima/m54081_181221_163846.bc1015_BAK8B_OA--bc1015_BAK8B_OA.bam",
        "resources/mock9/reads/lima/m54081_181221_163846.bc1016_BAK8B_OA--bc1016_BAK8B_OA.bam",
        "resources/mock9/reads/lima/m54081_181221_163846.bc1018_BAK8B_OA--bc1018_BAK8B_OA.bam",
        "resources/mock9/reads/lima/m54081_181221_163846.bc1019_BAK8B_OA--bc1019_BAK8B_OA.bam",
        "resources/mock9/reads/lima/m54081_181221_163846.bc1022_BAK8B_OA--bc1022_BAK8B_OA.bam"
    output:
        "resources/mock9/reads/lima/m54081_181221_163846.bc1001_BAK8A_OA--bc1001_BAK8A_OA.fq.gz",
        "resources/mock9/reads/lima/m54081_181221_163846.bc1002_BAK8A_OA--bc1002_BAK8A_OA.fq.gz",
        "resources/mock9/reads/lima/m54081_181221_163846.bc1009_BAK8A_OA--bc1009_BAK8A_OA.fq.gz",
        "resources/mock9/reads/lima/m54081_181221_163846.bc1010_BAK8A_OA--bc1010_BAK8A_OA.fq.gz",
        "resources/mock9/reads/lima/m54081_181221_163846.bc1012_BAK8A_OA--bc1012_BAK8A_OA.fq.gz",
        "resources/mock9/reads/lima/m54081_181221_163846.bc1015_BAK8B_OA--bc1015_BAK8B_OA.fq.gz",
        "resources/mock9/reads/lima/m54081_181221_163846.bc1016_BAK8B_OA--bc1016_BAK8B_OA.fq.gz",
        "resources/mock9/reads/lima/m54081_181221_163846.bc1018_BAK8B_OA--bc1018_BAK8B_OA.fq.gz",
        "resources/mock9/reads/lima/m54081_181221_163846.bc1019_BAK8B_OA--bc1019_BAK8B_OA.fq.gz",
        "resources/mock9/reads/lima/m54081_181221_163846.bc1022_BAK8B_OA--bc1022_BAK8B_OA.fq.gz"
    conda: 'envs/samtools.yaml'
    threads: 4
    shell:
        """
        for f in {input}; do
            fbase=${{f##*/}}
            fname=${{fbase%.bam}}
            samtools fastq --threads {threads} ${{f}} | gzip -dc >resources/mock9/reads/${{fname}}.fq.gz
        done
        """

#####################
### MOCK3 DATASET ###
#####################


####################
### NWC2 DATASET ###
####################

rule nwc2_pacbio_reads:
    input:  ancient("config/nwc2_pacbio.yaml")
    output: "resources/nwc2_pacbio/source/SRR7585901.fastq.gz"
    conda:  "envs/ncbi.yaml"
    shell:  "fasterq-dump -O resources/nwc2_pacbio/source -o SRR7585901.fastq --min-read-len 5000 SRR7585901 && gzip resources/nwc2_pacbio/source/SRR7585901.fastq"

#rule nwc2_pacbio_references:
#    output:
#        "resources/nwc2_pacbio/references/sthermophilus_nwc_2_1.fa",
#        "resources/nwc2_pacbio/references/ldelbrueckii_nwc_2_2.fa",
#        "resources/nwc2_pacbio/references/lhelveticus_nwc_2_3.fa",
#        "resources/nwc2_pacbio/references/lhelveticus_nwc_2_4.fa",
#        "resources/nwc2_pacbio/references.csv"
#    conda: "envs/ncbi.yaml"
#    shell:
#        """
#        (esearch -db nucleotide -query "CP031021.1" | efetch -format fasta >{output[0]}) && sleep 1 && \
#        (esearch -db nucleotide -query "CP031023.1" | efetch -format fasta >{output[1]}) && sleep 1 && \
#        (esearch -db nucleotide -query "CP031016.1" | efetch -format fasta >{output[2]}) && sleep 1 && \
#        (esearch -db nucleotide -query "CP031018.1" | efetch -format fasta >{output[3]}) && \
#        (for fpath in {output}; do fbase=${{fpath##*/}}; printf "%s,%s\n" "${{fbase%.fa}}" "$(readlink -f "${{fpath}}")"; done >{output[4]})
#        """

rule nwc2_ont_reads:
    input:  ancient("config/nwc2_ont.yaml")
    output: "resources/nwc2_ont/source/SRR7585900.fastq.gz"
    conda:  "envs/ncbi.yaml"
    shell:  "fasterq-dump -O resources/nwc2_ont/source -o SRR7585900.fastq --min-read-len 10000 SRR7585900 && gzip resources/nwc2_ont/source/SRR7585900.fastq"


###################
### HSM DATASET ###
###################



