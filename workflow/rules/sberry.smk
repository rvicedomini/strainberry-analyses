
rule strainberry:
    input:
        fasta="resources/{sample}/assembly/{assembly}.fa",
        faidx="resources/{sample}/assembly/{assembly}.fa.fai",
        bam="resources/{sample}/alignment/{assembly}.bam",
        bamidx="resources/{sample}/alignment/{assembly}.bam.bai"
    output:
        "results/{sample}/sberry_{assembly}_n{nstrains}/assembly.contigs.fa",
        "results/{sample}/sberry_{assembly}_n{nstrains}/assembly.scaffolds.fa"
    params:
        nstrains=config['nstrains']
    log:
        'logs/{sample}/sberry_{assembly}_n{nstrains}.log'
    conda:
        "../envs/sberry.yaml"
    threads: 14
    shell:
        "strainberry -r {input.fasta} -b {input.bam} -o sberry_{wildcards.assembly}_n{wildcards.nstrains} -n {wildcards.nstrains} -t {threads} &>{log}"
