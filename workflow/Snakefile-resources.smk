# The main entry point of your workflow.
# After configuring, running snakemake -n in a clone of this repository should successfully execute a dry-run of the workflow.
#report: "report/workflow.rst"

localrules: 
    dl_references,
    mock9_dl_data, mock9_dl_barcodes, 
    nwc2_dl_pacbio_reads, nwc2_dl_ont_reads,
    hsm_dl_reads, hsm_dl_assembly, hsm_dl_references

# Common rules, such as fasta/bam indexing, etc.
#include: "rules/common.smk"

# Allow users to fix the underlying OS via singularity.
#singularity: "docker://continuumio/miniconda3" # Allow users to fix the underlying OS via singularity.


rule all:
    input:
        # Mock datasets
        "resources/mock9/reads/m54081_181221_163846.bc1001_BAK8A_OA--bc1001_BAK8A_OA.fq.gz",
        "resources/mock9/reads/m54081_181221_163846.bc1002_BAK8A_OA--bc1002_BAK8A_OA.fq.gz",
        "resources/mock9/reads/m54081_181221_163846.bc1009_BAK8A_OA--bc1009_BAK8A_OA.fq.gz",
        "resources/mock9/reads/m54081_181221_163846.bc1010_BAK8A_OA--bc1010_BAK8A_OA.fq.gz",
        "resources/mock9/reads/m54081_181221_163846.bc1012_BAK8A_OA--bc1012_BAK8A_OA.fq.gz",
        "resources/mock9/reads/m54081_181221_163846.bc1015_BAK8B_OA--bc1015_BAK8B_OA.fq.gz",
        "resources/mock9/reads/m54081_181221_163846.bc1016_BAK8B_OA--bc1016_BAK8B_OA.fq.gz",
        "resources/mock9/reads/m54081_181221_163846.bc1018_BAK8B_OA--bc1018_BAK8B_OA.fq.gz",
        "resources/mock9/reads/m54081_181221_163846.bc1019_BAK8B_OA--bc1019_BAK8B_OA.fq.gz",
        "resources/mock9/reads/m54081_181221_163846.bc1022_BAK8B_OA--bc1022_BAK8B_OA.fq.gz",
        "resources/mock9/references.csv",
        "resources/mock3/reads/m54081_181221_163846.bc1002_BAK8A_OA--bc1002_BAK8A_OA.fq.gz",
        "resources/mock3/reads/m54081_181221_163846.bc1010_BAK8A_OA--bc1010_BAK8A_OA.fq.gz",
        "resources/mock3/reads/m54081_181221_163846.bc1019_BAK8B_OA--bc1019_BAK8B_OA.fq.gz",
        "resources/mock3/references.csv",
        # NWC2 - PacBio
        "resources/nwc2_pacbio/reads/SRR7585901.fastq.gz",
        "resources/nwc2_pacbio/references.csv",
        # NWC2 - Nanopore
        "resources/nwc2_ont/reads/SRR7585900.fastq.gz",
        "resources/nwc2_ont/references.csv",
        # Human Stool Microbiome
        "resources/hsm/reads/SRR8427258.fastq.gz",
        "resources/hsm/assembly/GCA_011075405.1.fasta",
        "resources/hsm/references.csv",


#########################
### REFERENCE GENOMES ###
#########################

mock9_reflist = [
    "NZ_CP034551.1", # B. cereus ATCC 14579
    "NC_000913.3",   # E. coli str. K-12 substr. MG1655
    "NC_017664.1",   # E. coli W ATCC 9637
    "NZ_CP006659.2", # K. pneumoniae ATCC BAA-2146
    "NZ_CP011398.2", # L. monocytogenes CFSAN008100
    "NC_008767.1",   # N. meningitidis FAM18
    "NZ_CP009361.1", # S. aureus ATCC 25923
    "NZ_CP041010.1", # S. aureus FDAARGOS 766
    "NZ_CP023645.1", # S. sonnei CFSAN030807
]

mock3_reflist = mock9_reflist[0:3]

nwc2_reflist = [
    "CP031021.1", # S. thermophilus NWC_2_1
    "CP031023.1", # L. delbrueckii NWC_2_2
    "CP031016.1", # L. helveticus NWC_2_3
    "CP031018.1", # L. helveticus NWC_2_4
]

all_reflist = mock9_reflist + nwc2_reflist

rule dl_references:
    output: 
        expand("resources/references/{acc}.fasta", acc=all_reflist)
    conda: 
        'envs/ncbi.yaml'
    shell: 
        "mkdir -p resources/references && efetch -format fasta -db nucleotide -id \"" + ",".join(all_reflist) + "\""
        "| awk 'BEGIN{{sid=\"undefined\"}} /^>/{{sid=substr($1,2)}} {{ print >\"resources/references/\"sid\".fasta\" }}'"


#####################
### MOCK9 DATASET ###
#####################

rule mock9_dl_data:
    input: ancient("config/mock9.yaml")
    output: "resources/mock9/source/MM_10Plex_42MB_SEQUENCING_DATA_EXPRESS2.0.tar.gz"
    shell: 
        """
        curl --create-dirs -L -o resources/mock9/source/MM_10Plex_42MB_SEQUENCING_DATA_EXPRESS2.0.tar.gz \
                https://downloads.pacbcloud.com/public/dataset/microbial_multiplex_dataset_release_SMRT_Link_v6.0.0_with_Express_2.0/MM_10Plex_42MB_SEQUENCING_DATA_EXPRESS2.0.tar.gz
        """

rule mock9_dl_barcodes:
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
        "resources/mock9/reads/demultiplexed/m54081_181221_163846.bc1001_BAK8A_OA--bc1001_BAK8A_OA.bam",
        "resources/mock9/reads/demultiplexed/m54081_181221_163846.bc1002_BAK8A_OA--bc1002_BAK8A_OA.bam",
        "resources/mock9/reads/demultiplexed/m54081_181221_163846.bc1009_BAK8A_OA--bc1009_BAK8A_OA.bam",
        "resources/mock9/reads/demultiplexed/m54081_181221_163846.bc1010_BAK8A_OA--bc1010_BAK8A_OA.bam",
        "resources/mock9/reads/demultiplexed/m54081_181221_163846.bc1012_BAK8A_OA--bc1012_BAK8A_OA.bam",
        "resources/mock9/reads/demultiplexed/m54081_181221_163846.bc1015_BAK8B_OA--bc1015_BAK8B_OA.bam",
        "resources/mock9/reads/demultiplexed/m54081_181221_163846.bc1016_BAK8B_OA--bc1016_BAK8B_OA.bam",
        "resources/mock9/reads/demultiplexed/m54081_181221_163846.bc1018_BAK8B_OA--bc1018_BAK8B_OA.bam",
        "resources/mock9/reads/demultiplexed/m54081_181221_163846.bc1019_BAK8B_OA--bc1019_BAK8B_OA.bam",
        "resources/mock9/reads/demultiplexed/m54081_181221_163846.bc1022_BAK8B_OA--bc1022_BAK8B_OA.bam"
    conda: 'envs/pacbio.yaml'
    threads: 
        workflow.cores
    shell:
        """
        tar -xf resources/mock9/source/MM_10Plex_42MB_SEQUENCING_DATA_EXPRESS2.0.tar.gz -C resources/mock9/source \
                && mkdir -p resources/mock9/reads/demultiplexed \
                && lima --num-threads {threads} --same --split-named \
                    resources/mock9/source/MM_10Plex_42MB_SEQUENCING_DATA_EXPRESS2.0/m54081_181221_163846.subreads.bam \
                    resources/mock9/source/Sequel_16_Barcodes_v3.fasta \
                    resources/mock9/reads/demultiplexed/m54081_181221_163846.bam
        """

rule mock9_fastq:
    input:
        ancient("resources/mock9/reads/demultiplexed/m54081_181221_163846.bc1001_BAK8A_OA--bc1001_BAK8A_OA.bam"),
        ancient("resources/mock9/reads/demultiplexed/m54081_181221_163846.bc1002_BAK8A_OA--bc1002_BAK8A_OA.bam"),
        ancient("resources/mock9/reads/demultiplexed/m54081_181221_163846.bc1009_BAK8A_OA--bc1009_BAK8A_OA.bam"),
        ancient("resources/mock9/reads/demultiplexed/m54081_181221_163846.bc1010_BAK8A_OA--bc1010_BAK8A_OA.bam"),
        ancient("resources/mock9/reads/demultiplexed/m54081_181221_163846.bc1012_BAK8A_OA--bc1012_BAK8A_OA.bam"),
        ancient("resources/mock9/reads/demultiplexed/m54081_181221_163846.bc1015_BAK8B_OA--bc1015_BAK8B_OA.bam"),
        ancient("resources/mock9/reads/demultiplexed/m54081_181221_163846.bc1016_BAK8B_OA--bc1016_BAK8B_OA.bam"),
        ancient("resources/mock9/reads/demultiplexed/m54081_181221_163846.bc1018_BAK8B_OA--bc1018_BAK8B_OA.bam"),
        ancient("resources/mock9/reads/demultiplexed/m54081_181221_163846.bc1019_BAK8B_OA--bc1019_BAK8B_OA.bam"),
        ancient("resources/mock9/reads/demultiplexed/m54081_181221_163846.bc1022_BAK8B_OA--bc1022_BAK8B_OA.bam")
    output:
        "resources/mock9/reads/m54081_181221_163846.bc1001_BAK8A_OA--bc1001_BAK8A_OA.fq.gz",
        "resources/mock9/reads/m54081_181221_163846.bc1002_BAK8A_OA--bc1002_BAK8A_OA.fq.gz",
        "resources/mock9/reads/m54081_181221_163846.bc1009_BAK8A_OA--bc1009_BAK8A_OA.fq.gz",
        "resources/mock9/reads/m54081_181221_163846.bc1010_BAK8A_OA--bc1010_BAK8A_OA.fq.gz",
        "resources/mock9/reads/m54081_181221_163846.bc1012_BAK8A_OA--bc1012_BAK8A_OA.fq.gz",
        "resources/mock9/reads/m54081_181221_163846.bc1015_BAK8B_OA--bc1015_BAK8B_OA.fq.gz",
        "resources/mock9/reads/m54081_181221_163846.bc1016_BAK8B_OA--bc1016_BAK8B_OA.fq.gz",
        "resources/mock9/reads/m54081_181221_163846.bc1018_BAK8B_OA--bc1018_BAK8B_OA.fq.gz",
        "resources/mock9/reads/m54081_181221_163846.bc1019_BAK8B_OA--bc1019_BAK8B_OA.fq.gz",
        "resources/mock9/reads/m54081_181221_163846.bc1022_BAK8B_OA--bc1022_BAK8B_OA.fq.gz"
    conda: 'envs/samtools.yaml'
    threads: 4
    shell:
        """
        for f in {input}; do
            fbase=${{f##*/}}
            fname=${{fbase%.bam}}
            samtools fastq --threads {threads} ${{f}} | gzip -c >resources/mock9/reads/${{fname}}.fq.gz
        done
        """

rule mock9_reflist:
    input:
        expand("resources/references/{acc}.fasta", acc=mock9_reflist)
    output: 
        "resources/mock9/references.csv"
    shell:
        """
        mkdir -p resources/mock9
        : > {output}
        printf "%s,%s\\n" "Bcereus" "$(readlink -f resources/references/NZ_CP034551.1.fasta)" >>{output}
        printf "%s,%s\\n" "Ecoli-K12" "$(readlink -f resources/references/NC_000913.3.fasta)" >>{output}
        printf "%s,%s\\n" "Ecoli-W" "$(readlink -f resources/references/NC_017664.1.fasta)" >>{output}
        printf "%s,%s\\n" "Kpneumoniae" "$(readlink -f resources/references/NZ_CP006659.2.fasta)" >>{output}
        printf "%s,%s\\n" "Lmonocytogenes" "$(readlink -f resources/references/NZ_CP011398.2.fasta)" >>{output}
        printf "%s,%s\\n" "Nmeningitidis" "$(readlink -f resources/references/NC_008767.1.fasta)" >>{output}
        printf "%s,%s\\n" "Saureus-ATCC" "$(readlink -f resources/references/NZ_CP009361.1.fasta)" >>{output}
        printf "%s,%s\\n" "Saureus-FDAA" "$(readlink -f resources/references/NZ_CP041010.1.fasta)" >>{output}
        printf "%s,%s\\n" "Ssonnei" "$(readlink -f resources/references/NZ_CP023645.1.fasta)" >>{output}
        """


#####################
### MOCK3 DATASET ###
#####################

rule mock3_fastq:
    input:
        ancient("resources/mock9/reads/m54081_181221_163846.bc1002_BAK8A_OA--bc1002_BAK8A_OA.fq.gz"),
        ancient("resources/mock9/reads/m54081_181221_163846.bc1010_BAK8A_OA--bc1010_BAK8A_OA.fq.gz"),
        ancient("resources/mock9/reads/m54081_181221_163846.bc1019_BAK8B_OA--bc1019_BAK8B_OA.fq.gz"),
    output:
        "resources/mock3/reads/m54081_181221_163846.bc1002_BAK8A_OA--bc1002_BAK8A_OA.fq.gz",
        "resources/mock3/reads/m54081_181221_163846.bc1010_BAK8A_OA--bc1010_BAK8A_OA.fq.gz",
        "resources/mock3/reads/m54081_181221_163846.bc1019_BAK8B_OA--bc1019_BAK8B_OA.fq.gz",
    shell:
        """
        mkdir -p resources/mock3/reads \
            && ln -s "$(readlink -f {input[0]})" {output[0]} \
            && ln -s "$(readlink -f {input[1]})" {output[1]} \
            && ln -s "$(readlink -f {input[2]})" {output[2]}
        """

rule mock3_reflist:
    input:
        expand("resources/references/{acc}.fasta", acc=mock3_reflist)
    output: 
        "resources/mock3/references.csv"
    shell:
        """
        mkdir -p resources/mock3
        : > {output}
        printf "%s,%s\\n" "Bcereus" "$(readlink -f resources/references/NZ_CP034551.1.fasta)" >>{output}
        printf "%s,%s\\n" "Ecoli-K12" "$(readlink -f resources/references/NC_000913.3.fasta)" >>{output}
        printf "%s,%s\\n" "Ecoli-W" "$(readlink -f resources/references/NC_017664.1.fasta)" >>{output}
        """

####################
### NWC2 DATASET ###
####################


rule nwc2_dl_pacbio_reads:
    input:  ancient("config/nwc2_pacbio.yaml")
    output: "resources/nwc2_pacbio/reads/SRR7585901.fastq.gz"
    conda:  "envs/ncbi.yaml"
    shell:  "fasterq-dump -O resources/nwc2_pacbio/reads -o SRR7585901.fastq --min-read-len 5000 SRR7585901 && gzip resources/nwc2_pacbio/reads/SRR7585901.fastq"

rule nwc2_pacbio_reflist:
    input:
        expand("resources/references/{acc}.fasta", acc=nwc2_reflist)
    output: 
        "resources/nwc2_pacbio/references.csv"
    shell:
        """
        : > {output}
        printf "%s,%s\\n" "sthermophilus_nwc2_1" "$(readlink -f resources/references/CP031021.1.fasta)" >>{output}
        printf "%s,%s\\n" "ldelbrueckii_nwc2_2"  "$(readlink -f resources/references/CP031023.1.fasta)" >>{output}
        printf "%s,%s\\n" "lhelveticus_nwc2_3"   "$(readlink -f resources/references/CP031016.1.fasta)" >>{output}
        printf "%s,%s\\n" "lhelveticus_nwc2_4"   "$(readlink -f resources/references/CP031018.1.fasta)" >>{output}
        """

rule nwc2_dl_ont_reads:
    input:  ancient("config/nwc2_ont.yaml")
    output: "resources/nwc2_ont/reads/SRR7585900.fastq.gz"
    conda:  "envs/ncbi.yaml"
    shell:  "fasterq-dump -O resources/nwc2_ont/reads -o SRR7585900.fastq --min-read-len 10000 SRR7585900 && gzip resources/nwc2_ont/reads/SRR7585900.fastq"

rule nwc2_ont_reflist:
    input:
        expand("resources/references/{acc}.fasta", acc=nwc2_reflist)
    output: 
        "resources/nwc2_ont/references.csv"
    shell:
        """
        : > {output}
        printf "%s,%s\\n" "sthermophilus_nwc2_1" "$(readlink -f resources/references/CP031021.1.fasta)" >>{output}
        printf "%s,%s\\n" "ldelbrueckii_nwc2_2"  "$(readlink -f resources/references/CP031023.1.fasta)" >>{output}
        printf "%s,%s\\n" "lhelveticus_nwc2_3"   "$(readlink -f resources/references/CP031016.1.fasta)" >>{output}
        printf "%s,%s\\n" "lhelveticus_nwc2_4"   "$(readlink -f resources/references/CP031018.1.fasta)" >>{output}
        """

###################
### HSM DATASET ###
###################

rule hsm_dl_reads:
    input:  ancient("config/hsm.yaml")
    output: "resources/hsm/reads/SRR8427258.fastq.gz"
    conda:  "envs/ncbi.yaml"
    shell:  "fasterq-dump -O resources/hsm/reads -o SRR8427258.fastq SRR8427258 && gzip resources/hsm/reads/SRR8427258.fastq"

rule hsm_dl_assembly:
    input:  ancient('config/hsm.yaml')
    output: 'resources/hsm/assembly/GCA_011075405.1.fasta'
    conda:  'envs/ncbi.yaml'
    shell:  
        '''
        mkdir -p resources/hsm/assembly
        esearch -db assembly -query "GCA_011075405.1" \
                | efetch -format docsum \
                | xtract -pattern DocumentSummary -element FtpPath_GenBank \
                | xargs -n 1 bash -c 'curl -s -L "${{0}}/${{0##*/}}_genomic.fna.gz" | gzip -dc >{output}'
        '''

rule hsm_dl_references:
    input: 
        ancient("config/hsm.yaml")
    output: 
        "resources/references/NZ_AEDR00000000.1.fasta",
        "resources/references/NZ_AEDS00000000.1.fasta",
    conda: 
        'envs/ncbi.yaml'
    shell:
        """
        mkdir -p resources/references
        esearch -db nucleotide -query "NZ_AEDR01000001.1:NZ_AEDR01000063.1[PACC]" | efetch -format fasta > resources/references/NZ_AEDR00000000.1.fasta
        sleep 2 && esearch -db nucleotide -query "NZ_AEDS01000001.1:NZ_AEDS01000070.1[PACC]" | efetch -format fasta > resources/references/NZ_AEDS00000000.1.fasta
        """

rule hsm_reflist:
    input:
        "resources/references/NZ_AEDR00000000.1.fasta",
        "resources/references/NZ_AEDS00000000.1.fasta",
    output: 
        "resources/hsm/references.csv"
    shell:
        """
        mkdir -p resources/hsm
        : > {output}
        printf "%s,%s\\n" "Vatypica-ACS-049-V-Sch6" "$(readlink -f resources/references/NZ_AEDR00000000.1.fasta)" >>{output}
        printf "%s,%s\\n" "Vatypica-ACS-134-V-Col7a" "$(readlink -f resources/references/NZ_AEDS00000000.1.fasta)" >>{output}
        """
