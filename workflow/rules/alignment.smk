
rule align_and_downsample:
    input:
        assembly='results/{sample}/assemblies/{assembly}.fa',
        reads=expand('{basepath}/{filename}', basepath=config['read_basepath'], filename=config['read_filenames'])
    output:
        bam='results/{sample}/alignments/{assembly}.bam'
    log:
        'logs/{sample}/{assembly}_minimap2.log'
    benchmark:
        'benchmarks/{sample}/{assembly}.alignment.benchmark.txt'
    conda:
        '../envs/alignment.yaml'
    threads:
        workflow.cores
    params:
        preset = 'map-ont' if config['technology'] == 'nanopore' else 'map-pb',
        subsample = f'-s {config["subsample_frac"]}' if config['subsample'] else '',
        bamthreads = min(4,workflow.cores)
    shell:
        """
        mkdir -p results/{wildcards.sample}/alignments
        minimap2 -t {threads} -2 -ax {params.preset} {input.assembly} {input.reads} 2>{log} \
          | samtools view -bhu - \
          | samtools sort --threads {threads} -u \
          | samtools view {params.subsample} --threads {params.bamthreads} -bh -o {output.bam} -
        """

