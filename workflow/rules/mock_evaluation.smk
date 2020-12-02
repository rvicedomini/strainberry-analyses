
rule assembly_stats:
    input:
        fasta="results/{sample}/assemblies/{assembly}.fa",
        references=config['ref_csv'],
    output:
        directory('results/{sample}/evaluation/{assembly}'),
        "results/{sample}/evaluation/{assembly}.report.tsv",
    log:
        'logs/{sample}/{assembly}_stats.log'
    conda:
        "../envs/evaluation.yaml"
    threads:
        workflow.cores
    shell:
        """
        python3 workflow/scripts/assembly_stats.py -f {input.fasta} -r {input.references} -o {output[0]} -t {threads} &>{log}
        cp {output[0]}/report.tsv {output[1]}
        """

