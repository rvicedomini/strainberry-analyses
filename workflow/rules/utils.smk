# FASTA index generation
rule samtools_faidx:
    input:  '{fname}.{faext}'
    output: '{fname}.{faext}.fai'
    wildcard_constraints: faext="fa|fasta|fna"
    conda:  '../envs/alignment.yaml'
    shell:  'samtools faidx {input}'

# BAM index generation
rule samtools_index:
    input:   '{fname}.bam'
    output:  '{fname}.bam.bai'
    conda:   '../envs/alignment.yaml'
    threads: 4
    shell:   'samtools index -@ {threads} {input}'

