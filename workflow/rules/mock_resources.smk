##############################
### MOCK REFERENCE GENOMES ###
##############################

mock_reflist = [
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

rule mock_dl_references:
    output: 
        expand("resources/mock/references/{acc}.fasta", acc=mock_reflist)
    conda: 
        '../envs/ncbi.yaml'
    shell: 
        "mkdir -p resources/mock/references && efetch -format fasta -db nucleotide -id \"" + ",".join(mock_reflist) + "\""
        "| awk 'BEGIN{{sid=\"undefined\"}} /^>/{{sid=substr($1,2)}} {{ print >\"resources/mock/references/\"sid\".fasta\" }}'"


#####################
### MOCK DATASET ###
#####################

rule mock_dl_data:
    output: 
        protected("resources/mock/source/MM_10Plex_42MB_SEQUENCING_DATA_EXPRESS2.0.tar.gz")
    shell: 
        """
        curl --create-dirs -L -o resources/mock/source/MM_10Plex_42MB_SEQUENCING_DATA_EXPRESS2.0.tar.gz \
          https://downloads.pacbcloud.com/public/dataset/microbial_multiplex_dataset_release_SMRT_Link_v6.0.0_with_Express_2.0/MM_10Plex_42MB_SEQUENCING_DATA_EXPRESS2.0.tar.gz
        """

rule mock_dl_barcodes:
    output: 
        protected("resources/mock/source/Sequel_16_Barcodes_v3.fasta")
    shell:  
        """
        curl --create-dirs -s -L -o resources/mock/source/Sequel_16_Barcodes_v3.zip https://www.pacb.com/wp-content/uploads/Sequel_16_Barcodes_v3.zip \
          && unzip -d resources/mock/source/ resources/mock/source/Sequel_16_Barcodes_v3.zip \
          && rm resources/mock/source/Sequel_16_Barcodes_v3.zip
        """

rule mock_demultiplex:
    input: 
        ancient("resources/mock/source/MM_10Plex_42MB_SEQUENCING_DATA_EXPRESS2.0.tar.gz"),
        ancient("resources/mock/source/Sequel_16_Barcodes_v3.fasta")
    output:
        "resources/mock/reads/demultiplexed/m54081_181221_163846.bc1001_BAK8A_OA--bc1001_BAK8A_OA.bam",
        "resources/mock/reads/demultiplexed/m54081_181221_163846.bc1002_BAK8A_OA--bc1002_BAK8A_OA.bam",
        "resources/mock/reads/demultiplexed/m54081_181221_163846.bc1009_BAK8A_OA--bc1009_BAK8A_OA.bam",
        "resources/mock/reads/demultiplexed/m54081_181221_163846.bc1010_BAK8A_OA--bc1010_BAK8A_OA.bam",
        "resources/mock/reads/demultiplexed/m54081_181221_163846.bc1012_BAK8A_OA--bc1012_BAK8A_OA.bam",
        "resources/mock/reads/demultiplexed/m54081_181221_163846.bc1015_BAK8B_OA--bc1015_BAK8B_OA.bam",
        "resources/mock/reads/demultiplexed/m54081_181221_163846.bc1016_BAK8B_OA--bc1016_BAK8B_OA.bam",
        "resources/mock/reads/demultiplexed/m54081_181221_163846.bc1018_BAK8B_OA--bc1018_BAK8B_OA.bam",
        "resources/mock/reads/demultiplexed/m54081_181221_163846.bc1019_BAK8B_OA--bc1019_BAK8B_OA.bam",
        "resources/mock/reads/demultiplexed/m54081_181221_163846.bc1022_BAK8B_OA--bc1022_BAK8B_OA.bam"
    conda: 
        '../envs/pacbio.yaml'
    threads: 
        workflow.cores
    shell:
        """
        tar -xf resources/mock/source/MM_10Plex_42MB_SEQUENCING_DATA_EXPRESS2.0.tar.gz -C resources/mock/source \
          && mkdir -p resources/mock/reads/demultiplexed \
          && lima --num-threads {threads} --split-bam-named --peek-guess \
               resources/mock/source/MM_10Plex_42MB_SEQUENCING_DATA_EXPRESS2.0/m54081_181221_163846.subreadset.xml \
               resources/mock/source/Sequel_16_Barcodes_v3.fasta \
               resources/mock/reads/demultiplexed/m54081_181221_163846.bam \
          && rm -r resources/mock/source/MM_10Plex_42MB_SEQUENCING_DATA_EXPRESS2.0
        """

rule mock_fastq:
    input:  
        ancient('resources/mock/reads/demultiplexed/{fname}.bam')
    output: 
        'resources/mock/reads/{fname}.fq.gz'
    conda: 
        '../envs/alignment.yaml'
    threads: 
        4
    shell: 
        'samtools fastq --threads {threads} {input} | gzip -c >{output}'


rule mock_reflist:
    input:
        expand("resources/mock/references/{acc}.fasta", acc=mock_reflist)
    output: 
        mock3="resources/mock/mock3.references.csv",
        mock9="resources/mock/mock9.references.csv",
    shell:
        """
        mkdir -p resources/mock && : > {output.mock3} && : > {output.mock9} \
          && (printf "%s,%s\\n" "Bcereus" "$(readlink -f resources/references/NZ_CP034551.1.fasta)" | tee --append {output.mock9} {output.mock3} >/dev/null) \
          && (printf "%s,%s\\n" "Ecoli-K12" "$(readlink -f resources/references/NC_000913.3.fasta)" | tee --append {output.mock9} {output.mock3} >/dev/null) \
          && (printf "%s,%s\\n" "Ecoli-W" "$(readlink -f resources/references/NC_017664.1.fasta)" | tee --append {output.mock9} {output.mock3} >/dev/null) \
          && (printf "%s,%s\\n" "Kpneumoniae" "$(readlink -f resources/references/NZ_CP006659.2.fasta)" | tee --append {output.mock9} >/dev/null) \
          && (printf "%s,%s\\n" "Lmonocytogenes" "$(readlink -f resources/references/NZ_CP011398.2.fasta)" | tee --append {output.mock9} >/dev/null) \
          && (printf "%s,%s\\n" "Nmeningitidis" "$(readlink -f resources/references/NC_008767.1.fasta)" | tee --append {output.mock9} >/dev/null) \
          && (printf "%s,%s\\n" "Saureus-ATCC" "$(readlink -f resources/references/NZ_CP009361.1.fasta)" | tee --append {output.mock9} >/dev/null) \
          && (printf "%s,%s\\n" "Saureus-FDAA" "$(readlink -f resources/references/NZ_CP041010.1.fasta)" | tee --append {output.mock9} >/dev/null) \
          && (printf "%s,%s\\n" "Ssonnei" "$(readlink -f resources/references/NZ_CP023645.1.fasta)" | tee --append {output.mock9} >/dev/null)
        """


