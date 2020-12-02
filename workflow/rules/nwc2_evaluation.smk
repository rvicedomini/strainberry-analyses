
rule nwc2_assembly_stats:
    input:
        fasta="results/{sample}/assemblies/{assembly}.fa",
        references=config['ref_csv'],
    output:
        "results/{sample}/evaluation/{assembly}.report.tsv",
        'results/{sample}/evaluation/{assembly}/assembly.Sthermophilus_NWC_2_1.fa',
        'results/{sample}/evaluation/{assembly}/assembly.Ldelbrueckii_NWC_2_2.fa',
        'results/{sample}/evaluation/{assembly}/assembly.Lhelveticus_NWC_2_3.fa',
        'results/{sample}/evaluation/{assembly}/assembly.Lhelveticus_NWC_2_4.fa',
    log:
        'logs/{sample}/{assembly}_stats.log'
    conda:
        "../envs/evaluation.yaml"
    threads:
        workflow.cores
    shell:
        """
        python3 workflow/scripts/assembly_stats.py -f {input.fasta} -r {input.references} -o results/{wildcards.sample}/evaluation/{wildcards.assembly} -t {threads} &>{log}
        cp results/{wildcards.sample}/evaluation/{wildcards.assembly}/report.tsv {output[0]}
        """

rule nwc2_checkm_bacteria:
    output: 'results/{sample}/evaluation/checkm/bacteria.ms'
    conda:  '../envs/evaluation.yaml'
    shell:  'mkdir -p results/{wildcards.sample}/evaluation/checkm && checkm taxon_set domain Bacteria {output[0]}'


rule nwc2_checkm_stats:
    input:
        asmdir=directory('results/{sample}/evaluation/{assembly}'),
        asmrep='results/{sample}/evaluation/{assembly}.report.tsv',
        bacms='results/{sample}/evaluation/checkm/bacteria.ms',
    output:
        checkm=directory('results/{sample}/evaluation/checkm/{assembly}'),
        tsv='results/{sample}/evaluation/{assembly}.checkm.tsv',
    log:
        'logs/{sample}/{assembly}_checkm.log'
    conda:
        "../envs/evaluation.yaml"
    threads: 
        8
    shell:
        """
        mkdir -p {output.checkm} \
          && checkm analyze -x fa -t {threads} {input.bacms} {input.asmdir} {output.checkm} 2>{log} \
          && checkm qa --tab_table -f {output.tsv} {input.bacms} {output.checkm} 2>>{log}
        """

