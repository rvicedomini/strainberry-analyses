
rule strainberry:
    input:
        fasta="results/{sample}/assemblies/{assembly}.fa",
        faidx="results/{sample}/assemblies/{assembly}.fa.fai",
        bam="results/{sample}/alignments/{assembly}.bam",
        bamidx="results/{sample}/alignments/{assembly}.bam.bai",
    output:
        sberry_ctg="results/{{sample}}/assemblies/sberry_{{assembly}}_n{0}_ctg.fa".format(config['nstrains']),
        sberry_scf="results/{{sample}}/assemblies/sberry_{{assembly}}_n{0}_scf.fa".format(config['nstrains']),
    params:
        nstrains=config['nstrains'],
        tech = '--nanopore' if config['technology'] == 'nanopore' else '',
    log:
        'logs/{sample}/sberry_{assembly}.log'
    conda:
        "../envs/strainberry.yaml"
    threads: 14
    shell:
        """
        strainberry -r {input.fasta} -b {input.bam} -o results/{wildcards.sample}/assemblies/sberry_{wildcards.assembly}_n{params.nstrains} -n {params.nstrains} -t {threads} {params.tech} &>{log} \
          && ln -s $(readlink -f results/{wildcards.sample}/assemblies/sberry_{wildcards.assembly}_n{params.nstrains}/assembly.contigs.fa) {output.sberry_ctg} \
          && ln -s $(readlink -f results/{wildcards.sample}/assemblies/sberry_{wildcards.assembly}_n{params.nstrains}/assembly.scaffolds.fa) {output.sberry_scf}
        """

