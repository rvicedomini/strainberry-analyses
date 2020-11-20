# FASTA index generation
rule samtools_faidx:
    input:  '{fname}.f{asta}'
    output: '{fname}.f{asta}.fai'
    conda:  '../envs/samtools.yaml'
    shell:  'samtools faidx {input}'

# BAM index generation
rule samtools_index:
    input:   '{fname}.bam'
    output:  '{fname}.bam.bai'
    conda:   '../envs/samtools.yaml'
    threads: 4
    shell:   'samtools index -@ {threads} {input}'

