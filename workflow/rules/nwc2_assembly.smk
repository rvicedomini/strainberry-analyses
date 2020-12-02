
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
        mkdir -p results/{wildcards.sample}/assemblies \
          && flye --meta --min-overlap 3000 --iterations 4 --genome-size {params.mgsize} --out-dir results/{wildcards.sample}/assemblies/flye --threads {threads} {params.readtype} {input.reads} &>{log} \
          && cp results/{wildcards.sample}/assemblies/flye/assembly.fasta {output[0]}
        """

rule canu_assembly:
    input: 
        reads=expand('{basepath}/{filename}', basepath=config['read_basepath'], filename=config['read_filenames'])
    output:
        protected('results/{sample}/assemblies/canu.fa')
    log:    
        'logs/{sample}/canu.log'
    conda:   
        '../envs/assembly.yaml'
    threads:
        workflow.cores
    params:
        readtype = '-nanopore' if config['technology'] == 'nanopore' else '-pacbio',
        mgsize = config['mg_size'],
    shell:
        """
        mkdir -p results/{wildcards.sample}/assemblies
        canu -p assembly -d results/{wildcards.sample}/assemblies/canu/ genomeSize={params.mgsize} maxThreads={threads} useGrid=false {params.readtype} {input.reads} &>{log} \
          && cp results/{wildcards.sample}/assemblies/canu/assembly.contigs.fasta {output[0]}
        """

