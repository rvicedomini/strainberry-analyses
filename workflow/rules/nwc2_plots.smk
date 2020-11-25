
rule nwc2_assembly_align:
    input:
        asm='results/{sample}/assembly_eval/{assembly}/assembly.{refname}.fa',
        ref_csv=config['ref_csv'],
    output:
        bam='results/{sample}/alignments/{assembly}-{refname}.bam'
    log:
        'logs/{sample}/minimap2_{assembly}-{refname}.log'
    conda:
        '../envs/alignment.yaml'
    threads:
        workflow.cores
    shell:
        """
        refpath="$( awk -F, -v rid="{wildcards.refname}" '$1==rid{{printf("%s",$2);exit}}' "{input.ref_csv}" )"
        minimap2 -ax asm20 -t {threads} --cs "${{refpath}}" {input.asm} | samtools sort --threads {threads} -o {output.bam}
        """


def nwc2_maxcov(wildcards):
    if wildcards.refname == 'Ldelbrueckii_NWC_2_2':
        if wildcards.assembly == 'flye':
            return 2
        elif wildcards.assembly == 'canu':
            return 3
    return 1


rule nwc2_refcoverage_plot:
    input:
        bam1='results/{sample}/alignments/{assembly}-{refname}.bam',
        bai1='results/{sample}/alignments/{assembly}-{refname}.bam.bai',
        bam2='results/{sample}/alignments/sberry_{assembly}_n2_scf-{refname}.bam',
        bai2='results/{sample}/alignments/sberry_{assembly}_n2_scf-{refname}.bam.bai',
    output:
        bedgz1='results/{sample}/assembly_eval/mosdepth/{assembly}-{refname}.regions.bed.gz',
        bedgz2='results/{sample}/assembly_eval/mosdepth/sberry_{assembly}_n2_scf-{refname}.regions.bed.gz',
        svg='results/{sample}/assembly_eval/refcoverage_{refname}-{assembly}.svg'
    log:
        'logs/{sample}/refcoverage_{refname}-{assembly}.log'
    conda:
        '../envs/assembly_eval.yaml'
    threads:
        workflow.cores
    params:
        maxcov=nwc2_maxcov
    shell:
        """
        mkdir -p results/{wildcards.sample}/assembly_eval/mosdepth \
          && mosdepth -t {threads} --by 20000 -m -x results/{wildcards.sample}/assembly_eval/mosdepth/{wildcards.assembly}-{wildcards.refname} {input.bam1} >{log} \
          && mosdepth -t {threads} --by 20000 -m -x results/{wildcards.sample}/assembly_eval/mosdepth/sberry_{wildcards.assembly}_n2_scf-{wildcards.refname} {input.bam2} >{log} \
          && python3 workflow/scripts/assembly_refcoverage.py -m {params.maxcov} -b {output.bedgz1} {output.bedgz2} \
               -p results/{wildcards.sample}/assembly_eval/refcoverage_{wildcards.refname}-{wildcards.assembly} >>{log}
        """

