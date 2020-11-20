
rule assembly_stats:
    input:
        fasta="results/{sample}/assemblies/{assembly}.fa",
        references=config['ref_csv'],
    output:
        "results/{sample}/assembly_eval/{assembly}.report.tsv"
    log:
        'logs/{sample}/{assembly}_stats.log'
    conda:
        "../envs/assembly_eval.yaml"
    threads:
        workflow.cores
    shell:
        """
        python3 workflow/scripts/assembly_stats.py -f {input.fasta} -r {input.references} -o results/{sample}/assembly_eval/{wildcards.assembly} -t {threads} &>{log} \
        cp results/{sample}/assembly_eval/{wildcards.assembly}/report.tsv {output[0]}
        """

