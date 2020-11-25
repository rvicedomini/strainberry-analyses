
rule nwc2_assembly_stats:
    input:
        fasta="results/{sample}/assemblies/{assembly}.fa",
        references=config['ref_csv'],
    output:
        "results/{sample}/assembly_eval/{assembly}.report.tsv",
        'results/{sample}/assembly_eval/{assembly}/assembly.Sthermophilus_NWC_2_1.fa',
        'results/{sample}/assembly_eval/{assembly}/assembly.Ldelbrueckii_NWC_2_2.fa',
        'results/{sample}/assembly_eval/{assembly}/assembly.Lhelveticus_NWC_2_3.fa',
        'results/{sample}/assembly_eval/{assembly}/assembly.Lhelveticus_NWC_2_4.fa',
    log:
        'logs/{sample}/{assembly}_stats.log'
    conda:
        "../envs/assembly_eval.yaml"
    threads:
        workflow.cores
    shell:
        """
        python3 workflow/scripts/assembly_stats.py -f {input.fasta} -r {input.references} -o results/{sample}/assembly_eval/{wildcards.assembly} -t {threads} &>{log}
        cp results/{sample}/assembly_eval/{wildcards.assembly}/report.tsv {output[0]}
        """

rule nwc2_checkm_stats:
    input:
        'results/{sample}/assembly_eval/{assembly}.report.tsv',
    output:
        'results/{sample}/assembly_eval/{assembly}.checkm.tsv'
    log:
        'logs/{sample}/{assembly}_checkm.log'
    conda:
        "../envs/assembly_eval.yaml"
    threads:
        workflow.cores
    shell:
        """
        mkdir -p results/{sample}/assembly_eval/checkm_{assembly} \
          && checkm taxon_set domain Bacteria \
               results/{sample}/assembly_eval/checkm_{assembly}/bacteria.ms 2>{log} \
          && checkm analyze -x fa -t {threads} \
               results/{sample}/assembly_eval/checkm_{assembly}/bacteria.ms \
               results/{sample}/assembly_eval/{assembly}/ \
               results/{sample}/assembly_eval/checkm_{assembly} 2>>{log} \
          && checkm qa --tab_table -f {output[0]} \
               results/{sample}/assembly_eval/checkm_{assembly}/bacteria.ms \
               results/{sample}/assembly_eval/checkm_{assembly} 2>>{log}
        """

