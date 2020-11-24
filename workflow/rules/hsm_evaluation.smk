
rule hsm_kraken2:
    input:
        fasta='results/{sample}/assemblies/{assembly}.fa',
    output:
        kreport='results/{sample}/kraken2/{assembly}.kreport',
        kraken2='results/{sample}/kraken2/{assembly}.kraken2',
    conda:
        '../envs/assembly_eval.yaml',
    threads:
        workflow.cores
    params:
        db=config['kraken2_db']
    shell:
        """
        kraken2 --db {params.db} --threads {threads} --report {output.kreport} --output {output.kraken2} --use-names {input.fasta}
        """

rule hsm_metabat2:
    input:
        fasta='results/{sample}/assemblies/{assembly}.fa',
        faidx='results/{sample}/assemblies/{assembly}.fa.fai',
        bam='results/{sample}/alignments/{assembly}.bam',
        bamidx='results/{sample}/alignments/{assembly}.bam.bai',
    output:
        depth='results/{sample}/binning/{assembly}.depth.txt',
        bins=directory('results/{sample}/binning/{assembly}'),
    log:
        'logs/{sample}/{assembly}_metabat2.log'
    conda:
        "../envs/assembly_eval.yaml"
    threads:
        workflow.cores
    shell:
        """
        mkdir -p {output.bins} \
          && jgi_summarize_bam_contig_depths --minContigLength 1000 --minContigDepth 1 --percentIdentity 50 --outputDepth {output.depth} {input.bam} 2>{log} \
          && metabat2 -t {threads} -i {input.fasta} -a {output.depth} -o {output.bins}/bin --unbinned --saveCls --seed 1 2>>{log}
        """

rule hsm_checkm_eval:
    input:
        depth='results/{sample}/binning/{assembly}.depth.txt',
        bins=directory('results/{sample}/binning/{assembly}'),
    output:
        checkm='results/{sample}/assembly_eval/{assembly}.checkm.tsv'
    log:
        'logs/{sample}/{assembly}_checkm.log'
    conda:
        "../envs/assembly_eval.yaml"
    threads:
        workflow.cores
    shell:
        """
        checkm lineage_wf --tab_table -f {output.checkm} -t {threads} --pplacer_threads {threads} \
          -x fa {input.bins} results/{wildcards.sample}/binning/checkm_{wildcards.assembly}/
        """

