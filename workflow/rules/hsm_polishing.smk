# the following should not be necessary
#ruleorder: racon_init > racon_consensus

rule racon_paf:
    input:
        reads=expand('{basepath}/{filename}', basepath=config['read_basepath'], filename=config['read_filenames']),
        assembly='results/{sample}/assembly/racon/{assembly}.racon_{iter}.fa',
    output:
        pafgz='results/{sample}/assembly/racon/{assembly}.racon_{iter}.paf.gz'
    conda:
        '../envs/polishing.yaml'
    threads:
        workflow.cores
    params:
        preset = 'map-ont' if config['technology'] == 'nanopore' else 'map-pb',
    shell:
        """
        minimap2 -t {threads} -x {params.preset} {input.assembly} {input.reads} | gzip -c >{output.pafgz}
        """

rule racon_init:
    input:  'results/{sample}/assembly/{assembly}.fa',
    output: 'results/{sample}/assembly/racon/{assembly}.racon_0.fa',
    shell:  'ln -s "$(readlink -f {input})" {output}'

def prev_racon_output(wildcards):
    iter=int(wildcards.iter)
    if iter < 0:
        raise ValueError("racon iteration value should be >= 1")
    pafgz=f'results/{wildcards.sample}/assembly/racon/{wildcards.assembly}.racon_{iter-1}.paf.gz'
    assembly=f'results/{wildcards.sample}/assembly/racon/{wildcards.assembly}.racon_{iter-1}.fa'
    return pafgz,assembly

rule racon_consensus:
    input:
        expand('{basepath}/{filename}', basepath=config['read_basepath'], filename=config['read_filenames']),
        prev_racon_output
    output:
        'results/{sample}/assembly/racon/{assembly}.racon_{iter}.fa',
    conda:
        '../envs/polishing.yaml'
    threads: 
        workflow.cores
    shell:
        """
        racon -m 8 -x -6 -g -8 -w 500 -t {threads} {input} > {output[0]}
        """

rule racon_final:
    input:  'results/{sample}/assembly/racon/{assembly}.racon_{iter}.fa',
    output: 'results/{sample}/assembly/{assembly}.racon_{iter}.fa',
    shell:  'ln -s "$(readlink -f {input})" {output}'


checkpoint medaka_init:
    input:
        fasta='results/{sample}/assembly/{assembly}.racon_4.fa',
        faidx='results/{sample}/assembly/{assembly}.racon_4.fa.fai',
    output:
        ids=directory('results/{sample}/assembly/medaka/{assembly}/ids'),
        contigs=directory('results/{sample}/assembly/medaka/{assembly}/contigs'),
        results=directory('results/{sample}/assembly/medaka/{assembly}/results'),
        polished=directory('results/{sample}/assembly/medaka/{assembly}/polished'),
    conda:
        '../envs/polishing.yaml'
    shell:
        """
        mkdir -p {output.ids} {output.contigs} {output.results} {output.polished} \
          && (cut -f1 {input.faidx} | xargs -n 1 -I{{}} touch {output.ids}/{{}}) \
          && (cut -f1 {input.faidx} | xargs -n 1 -I{{}} samtools faidx {input.fasta} {{}} >{output.contigs}/{{}}.fa)
        """

rule medaka_consensus:
    input:
        'results/{sample}/assembly/medaka/{assembly}/ids/{contig}',
        'results/{sample}/alignments/{assembly}.racon_4.bam',
        'results/{sample}/alignments/{assembly}.racon_4.bam.bai',
    output:
        'results/{sample}/assembly/medaka/{assembly}/results/{contig}.hdf'
    conda:
        '../envs/polishing.yaml'
    shell:
        """
        medaka consensus {input[1]} {output} --model r941_flip213 --threads 1 --regions {wildcards.contig}:0-
        """

rule medaka_stitch:
    input:  'results/{sample}/assembly/2.polish/medaka/results/{contig}.hdf', #'results/{sample}/assembly/medaka/{assembly}/contigs{range}.fa',
    output: 'results/{sample}/assembly/medaka/{assembly}/polished/{contig}.fa'
    shell:  'medaka stitch {input[0]} {output[0]}'
    #shell:  'medaka stitch {input[0]} {input[1]} {output}'

def aggregate_medaka_results(wildcards):
    checkpoint_output = checkpoints.medaka_init.get(**wildcards).output[0]
    return expand('results/{sample}/assembly/medaka/{assembly}/polished/{contig}.fa',
            sample=wildcards.sample,
            assembly=wildcards.assembly,
            contig=glob_wildcards(os.path.join(checkpoint_output,'{contig}')).contig)

rule medaka_aggregate:
    input:
        aggregate_medaka_subsetruns
    output:
        'results/{sample}/assembly/{assembly}.medaka.fa'
    shell:
        """
        find results/{wildcards.sample}/assembly/medaka/{wildcards.assembly}/polished -maxdepth 1 -type f -name '*.fa' -print0 \
          | xargs -0 cat | awk '/^>/{ split($1,a,":"); printf("%s-%d\n",a[1],++x[a[1]]); next } {print}' \
          >{output[0]}
        """

