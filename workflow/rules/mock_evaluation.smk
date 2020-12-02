
rule assembly_stats:
    input:
        fasta="results/{sample}/assemblies/{assembly}.fa",
        references=config['ref_csv'],
    output:
        "results/{sample}/evaluation/{assembly}.report.tsv"
    log:
        'logs/{sample}/{assembly}_stats.log'
    conda:
        "../envs/evaluation.yaml"
    threads:
        workflow.cores
    shell:
        """
        python3 workflow/scripts/assembly_stats.py -f {input.fasta} -r {input.references} -o results/{sample}/evaluation/{wildcards.assembly} -t {threads} &>{log}
        cp results/{sample}/evaluation/{wildcards.assembly}/report.tsv {output[0]}
        """

