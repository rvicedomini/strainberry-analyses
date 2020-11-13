
rule assembly_eval:
    input:
        fasta="results/{sample}/assemblies/{assembly}.fa",
        references="resources/{sample}/references.csv",
    output:
        "results/{sample}/assembly_eval/{assembly}.report.tsv"
    log:
        'logs/{sample}/{assembly}_eval.log'
    conda:
        "../envs/assembly_eval.yaml"
    threads: 14
    shell:
        """
        python3 workflow/scripts/assembly_eval.py -f {input.fasta} -r {input.references} -o results/{sample}/assembly_eval/{wildcards.assembly} -t 14 &>{log}
        cp results/{sample}/assembly_eval/{wildcards.assembly}/report.tsv {output[0]}
        """
