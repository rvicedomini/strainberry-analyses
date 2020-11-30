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
    output: 
        "resources/hsm/references/NZ_AEDR00000000.1.fasta",
        "resources/hsm/references/NZ_AEDS00000000.1.fasta",
    conda: 
        'envs/ncbi.yaml'
    shell:
        """
        mkdir -p resources/hsm/references \
          && (esearch -db nucleotide -query "NZ_AEDR01000001.1:NZ_AEDR01000063.1[PACC]" | efetch -format fasta > resources/references/NZ_AEDR00000000.1.fasta) \
          && (esearch -db nucleotide -query "NZ_AEDS01000001.1:NZ_AEDS01000070.1[PACC]" | efetch -format fasta > resources/references/NZ_AEDS00000000.1.fasta)
        """

rule hsm_reflist:
    input:
        "resources/hsm/references/NZ_AEDR00000000.1.fasta",
        "resources/hsm/references/NZ_AEDS00000000.1.fasta",
    output: 
        "resources/hsm/vatypica.references.csv"
    shell:
        """
        : > {output[0]} \
          && (printf "%s,%s\\n" "Vatypica-ACS-049-V-Sch6" "$(readlink -f resources/hsm/references/NZ_AEDR00000000.1.fasta)"  | tee --append {output[0]} >/dev/null) \
          && (printf "%s,%s\\n" "Vatypica-ACS-134-V-Col7a" "$(readlink -f resources/hsm/references/NZ_AEDS00000000.1.fasta)" | tee --append {output[0]} >/dev/null)
        """

if config['download_kraken2_db']:
    rule hsm_dl_references:
        output:
            'resources/maxikraken2_1903_140GB/hash.k2d',
            'resources/maxikraken2_1903_140GB/opts.k2d',
            'resources/maxikraken2_1903_140GB/taxo.k2d',
        conda:
            'envs/ncbi.yaml'
        shell:
            """
            mkdir -p resources/maxikraken2_1903_140GB && cd resources/maxikraken2_1903_140GB \
              && wget --quiet -c https://refdb.s3.climb.ac.uk/maxikraken2_1903_140GB/hash.k2d \
              && wget --quiet https://refdb.s3.climb.ac.uk/maxikraken2_1903_140GB/opts.k2d \
              && wget --quiet https://refdb.s3.climb.ac.uk/maxikraken2_1903_140GB/taxo.k2d
            """

