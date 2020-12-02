rule racon_paf:
    input:
        reads=expand('{basepath}/{filename}', basepath=config['read_basepath'], filename=config['read_filenames']),
        assembly='results/{sample}/assemblies/racon/{assembly}.racon_{it}.fa',
    output:
        pafgz='results/{sample}/assemblies/racon/{assembly}.racon_{it}.paf.gz'
    conda:
        '../envs/alignment.yaml'
    threads:
        workflow.cores
    params:
        preset = 'map-ont' if config['technology'] == 'nanopore' else 'map-pb',
    shell:
        """
        minimap2 -t {threads} -x {params.preset} {input.assembly} {input.reads} | gzip -c >{output.pafgz}
        """

rule racon_init:
    input:  'results/{sample}/assemblies/{assembly}.fa',
    output: 'results/{sample}/assemblies/racon/{assembly}.racon_0.fa',
    shell:  'ln -s "$(readlink -f {input})" {output}'

def prev_racon_output(wildcards):
    it=int(wildcards.it)
    if it < 0:
        raise ValueError("racon iteration value should be >= 1")
    pafgz=f'results/{wildcards.sample}/assemblies/racon/{wildcards.assembly}.racon_{it-1}.paf.gz'
    assembly=f'results/{wildcards.sample}/assemblies/racon/{wildcards.assembly}.racon_{it-1}.fa'
    return (pafgz,assembly)

rule racon_consensus:
    input:
        expand('{basepath}/{filename}', basepath=config['read_basepath'], filename=config['read_filenames']),
        prev_racon_output,
    output:
        'results/{sample}/assemblies/racon/{assembly}.racon_{it}.fa',
    wildcard_constraints:
        it="\d+"
    conda:
        '../envs/polishing.yaml'
    threads: 
        workflow.cores
    shell:
        """
        racon -m 8 -x -6 -g -8 -w 500 -t {threads} {input} > {output[0]}
        """

rule racon_final:
    input:  'results/{sample}/assemblies/racon/{assembly}.racon_4.fa',
    output: 'results/{sample}/assemblies/racon/{assembly}.racon_4.renamed.fa',
    wildcard_constraints: assembly="[^/]+"
    shell:  'cat {input[0]} | sed "s/>sberry|/>sberry_/" >{output[0]}'

rule racon_alignment:
    input:
        fasta='results/{sample}/assemblies/racon/{assembly}.racon_4.renamed.fa',
        faidx='results/{sample}/assemblies/racon/{assembly}.racon_4.renamed.fa.fai',
        reads=expand('{basepath}/{filename}', basepath=config['read_basepath'], filename=config['read_filenames']),
    output:
        bam='results/{sample}/assemblies/racon/{assembly}.racon_4.renamed.bam',
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
        minimap2 -t {threads} -2 -ax {params.preset} {input.fasta} {input.reads} \
          | samtools sort --threads {threads} \
          | samtools view {params.subsample} --threads 4 -bSh -o {output.bam} -
        """

checkpoint medaka_init:
    input:
        fasta='results/{sample}/assemblies/racon/{assembly}.racon_4.renamed.fa',
        faidx='results/{sample}/assemblies/racon/{assembly}.racon_4.renamed.fa.fai',
    output:
        ids=directory('results/{sample}/assemblies/medaka/{assembly}/ids'),
        contigs=directory('results/{sample}/assemblies/medaka/{assembly}/contigs'),
        results=directory('results/{sample}/assemblies/medaka/{assembly}/results'),
        polished=directory('results/{sample}/assemblies/medaka/{assembly}/polished'),
    conda:
        '../envs/alignment.yaml'
    shell:
        """
        mkdir -p {output.ids} {output.contigs} {output.results} {output.polished} \
          && (cut -f1 {input.faidx} | xargs -n 1 -I{{}} touch '{output.ids}/{{}}') \
          && (cut -f1 {input.faidx} | xargs -n 1 -I{{}} bash -c 'samtools faidx {input.fasta} {{}} >{output.contigs}/{{}}.fa')
        """

rule medaka_consensus:
    input:
        'results/{sample}/assemblies/medaka/{assembly}/ids/{contig}',
        'results/{sample}/assemblies/racon/{assembly}.racon_4.renamed.bam',
        'results/{sample}/assemblies/racon/{assembly}.racon_4.renamed.bam.bai',
    output:
        'results/{sample}/assemblies/medaka/{assembly}/results/{contig}.hdf'
    conda:
        '../envs/polishing.yaml'
    shell:
        """
        medaka consensus {input[1]} '{output[0]}' --model r941_flip213 --threads 1 --regions '{wildcards.contig}:0-'
        """

rule medaka_stitch:
    input:  'results/{sample}/assemblies/medaka/{assembly}/results/{contig}.hdf',
    output: 'results/{sample}/assemblies/medaka/{assembly}/polished/{contig}.fa',
    conda:  '../envs/polishing.yaml'
    shell:  "medaka stitch '{input[0]}' '{output[0]}'"

def aggregate_medaka_results(wildcards):
    checkpoint_output = checkpoints.medaka_init.get(**wildcards).output.ids
    return expand('results/{sample}/assemblies/medaka/{assembly}/polished/{contig}.fa',
            sample=wildcards.sample,
            assembly=wildcards.assembly,
            contig=glob_wildcards(os.path.join(checkpoint_output,'{ctg}')).ctg)

rule medaka_aggregate:
    input:
        aggregate_medaka_results
    output:
        'results/{sample}/assemblies/{assembly}.medaka.fa'
    shell:
        """
        find results/{wildcards.sample}/assemblies/medaka/{wildcards.assembly}/polished -maxdepth 1 -type f -name '*.fa' -print0 \
          | xargs -0 cat | awk '/^>/{{ split($1,a,":"); printf("%s-%d\\n",a[1],++x[a[1]]); next }} {{print}}' \
          | sed 's#>sberry_#>sberry|#' >{output[0]}
        """

