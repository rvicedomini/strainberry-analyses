
rule strainberry:
    input:
        fasta="results/{sample}/assemblies/{assembly}.fa",
        faidx="results/{sample}/assemblies/{assembly}.fa.fai",
        bam="results/{sample}/alignments/{assembly}.bam",
        bamidx="results/{sample}/alignments/{assembly}.bam.bai",
    output:
        # "results/{sample}/assemblies/strain_separated/sberry_{assembly}_n{nstrains}/assembly.contigs.fa",
        # "results/{sample}/assemblies/strain_separated/sberry_{assembly}_n{nstrains}/assembly.scaffolds.fa",
        "results/{sample}/assemblies/sberry_{assembly}_n{nstrains}.contigs.fa",
        "results/{sample}/assemblies/sberry_{assembly}_n{nstrains}.scaffolds.fa",
    params:
        nstrains=config['nstrains']
    log:
        'logs/{sample}/sberry_{assembly}_n{nstrains}.log'
    conda:
        "../envs/sberry.yaml"
    threads: 14
    shell:
        """
        strainberry -r {input.fasta} -b {input.bam} -o results/{wildcards.sample}/assemblies/sberry_{wildcards.assembly}_n{wildcards.nstrains} -n {wildcards.nstrains} -t {threads} &>{log}
        cp results/{wildcards.sample}/assemblies/sberry_{wildcards.assembly}_n{wildcards.nstrains}/assembly.contigs.fa {output[0]}
        cp results/{wildcards.sample}/assemblies/sberry_{wildcards.assembly}_n{wildcards.nstrains}/assembly.scaffolds.fa {output[1]}
        """

