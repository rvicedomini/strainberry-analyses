
rule strainberry:
    input:
        fasta="results/{sample}/assemblies/{assembly}.fa",
        faidx="results/{sample}/assemblies/{assembly}.fa.fai",
        bam="results/{sample}/alignments/{assembly}.bam",
        bamidx="results/{sample}/alignments/{assembly}.bam.bai",
    output:
        sberry_asm="results/{sample}/assemblies/sberry_{assembly}.fa",
    params:
        tech = '--nanopore' if config['technology'] == 'nanopore' else '',
    log:
        'logs/{sample}/sberry_{assembly}.log'
    benchmark:
        'benchmarks/{sample}/sberry_{assembly}.benchmark.txt'
    conda:
        "../envs/strainberry.yaml"
    threads:
        workflow.cores
    shell:
        """
        strainberry -r {input.fasta} -b {input.bam} -o results/{wildcards.sample}/assemblies/sberry_{wildcards.assembly} -c {threads} {params.tech} &>{log} \
          && cp results/{wildcards.sample}/assemblies/sberry_{wildcards.assembly}/assembly.scaffolds.fa {output.sberry_asm}
        """

