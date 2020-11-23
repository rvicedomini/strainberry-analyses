#######################
### NWC2 REFERENCES ###
#######################

nwc2_reflist = [
    "CP031021.1", # S. thermophilus NWC_2_1
    "CP031023.1", # L. delbrueckii NWC_2_2
    "CP031016.1", # L. helveticus NWC_2_3
    "CP031018.1", # L. helveticus NWC_2_4
]

rule nwc2_dl_references:
    output: 
        expand("resources/nwc2/references/{acc}.fasta", acc=nwc2_reflist)
    conda: 
        'envs/ncbi.yaml'
    shell: 
        "mkdir -p resources/nwc2/references && efetch -format fasta -db nucleotide -id \"" + ",".join(nwc2_reflist) + "\""
        "| awk 'BEGIN{{sid=\"undefined\"}} /^>/{{sid=substr($1,2)}} {{ print >\"resources/nwc2/references/\"sid\".fasta\" }}'"

rule nwc2_reflist:
    input:
        expand("resources/nwc2/references/{acc}.fasta", acc=nwc2_reflist)
    output: 
        "resources/nwc2/nwc2.references.csv"
    shell:
        """
        : > {output} \
          && (printf "%s,%s\\n" "Sthermophilus_NWC_2_1" "$(readlink -f resources/references/CP031021.1.fasta)" | tee --append {output[0]} >/dev/null) \
          && (printf "%s,%s\\n" "Ldelbrueckii_NWC_2_2"  "$(readlink -f resources/references/CP031023.1.fasta)" | tee --append {output[0]} >/dev/null) \
          && (printf "%s,%s\\n" "Lhelveticus_NWC_2_3"   "$(readlink -f resources/references/CP031016.1.fasta)" | tee --append {output[0]} >/dev/null) \
          && (printf "%s,%s\\n" "Lhelveticus_NWC_2_4"   "$(readlink -f resources/references/CP031018.1.fasta)" | tee --append {output[0]} >/dev/null)
        """

####################
### NWC2 DATASET ###
####################

rule nwc2_pacbio_dl_reads:
    input:  ancient("config/nwc2_pacbio.yaml")
    output: "resources/nwc2_pacbio/reads/SRR7585901.fastq.gz"
    conda:  "envs/ncbi.yaml"
    shell:  "fasterq-dump -O resources/nwc2_pacbio/reads -o SRR7585901.fastq --min-read-len 5000 SRR7585901 && gzip resources/nwc2_pacbio/reads/SRR7585901.fastq"

rule nwc2_ont_dl_reads:
    input:  ancient("config/nwc2_ont.yaml")
    output: "resources/nwc2_ont/reads/SRR7585900.fastq.gz"
    conda:  "envs/ncbi.yaml"
    shell:  "fasterq-dump -O resources/nwc2_ont/reads -o SRR7585900.fastq --min-read-len 10000 SRR7585900 && gzip resources/nwc2_ont/reads/SRR7585900.fastq"


