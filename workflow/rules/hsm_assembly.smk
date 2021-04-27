
rule flye_assembly:
    input:
        reads=expand('{basepath}/{filename}', basepath=config['read_basepath'], filename=config['read_filenames'])
    output: 
        protected('results/{sample}/assemblies/flye.fa')
    log:    
        'logs/{sample}/flye.log'
    conda:   
        '../envs/assembly.yaml'
    threads:
        workflow.cores
    params:
        readtype = '--nano-raw' if config['technology'] == 'nanopore' else '--pacbio-raw',
        mgsize = config['mg_size']
    shell:
        """
        mkdir -p results/{wildcards.sample}/assemblies
        flye --meta --out-dir results/{wildcards.sample}/assemblies/flye --genome-size {params.mgsize} --threads {threads} {params.readtype} {input.reads} &>{log} \
          && cp results/{wildcards.sample}/assemblies/flye/assembly.fasta {output[0]}
        """

rule lathe_assembly:
    input:  'resources/hsm/assembly/GCA_011075405.1.fasta'
    output: 'results/{sample}/assemblies/lathe-p1.fa'
    shell:  'ln "$(readlink -f {input})" {output}'

