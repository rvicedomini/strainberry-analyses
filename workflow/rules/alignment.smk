
if config['read_has_qual']:

    rule align_and_downsample:
        input:
            assembly='results/{sample}/assemblies/{assembly}.fa',
            reads=expand('{basepath}/{filename}', basepath=config['read_basepath'], filename=config['read_filenames'])
        output:
            bam='results/{sample}/alignments/{assembly}.bam'
        log:
            'logs/{sample}/{assembly}_minimap2.log'
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
              | samtools sort --threads {threads} \
              | samtools view {params.subsample} --threads 4 -bSh -o {output.bam} -
            """

else:

    rule align_and_downsample:
        input:
            assembly='results/{sample}/assemblies/{assembly}.fa',
            reads=expand('{basepath}/{filename}', basepath=config['read_basepath'], filename=config['read_filenames'])
        output:
            bam='results/{sample}/alignments/{assembly}.bam'
        log:
            'logs/{sample}/{assembly}_minimap2.log'
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
              | samtools sort --threads {threads} \
              | samtools view {params.subsample} --threads {params.bamthreads} -bSh - \
              | python3 workflow/scripts/replace_base_qualities.py --threads {params.bamthreads} -b - -o {output.bam}
            """

