
checkpoint sberry_bins_init:
    input:
        lathe_bins='results/hsm/binning/lathe-p1',
        bin_stats='results/{sample}/evaluation/lathe-p1.bin_stats.tsv',
    output:
        best_bins=directory('results/{sample}/binning/lathe-p1_bestbins'),
        sberry_bins=directory('results/{sample}/binning/strainberry'),
        done='results/{sample}/binning/sberry_bestbins_init.done',
    shell:
        """
        mkdir -p {output.best_bins} {output.sberry_bins} \
          && (cut -f1 {input.bin_stats} | grep -Ev '^#' | sort | uniq | xargs -n 1 -I{{}} cp "{input.lathe_bins}/{{}}.fa" "{output.best_bins}") \
          && touch {output.done}
        """

rule sberry_checkm_eval:
    input:
        'results/{sample}/binning/sberry_bestbins_init.done',
        'results/{sample}/binning/strainberry',
    output:
        checkm='results/{sample}/binning/strainberry.checkm.tsv'
    log:
        'logs/{sample}/strainberry_checkm.log'
    conda:
        "../envs/evaluation.yaml"
    threads:
        workflow.cores
    shell:
        """
        checkm lineage_wf --tab_table -f {output.checkm} -t {threads} --pplacer_threads {threads} \
          -x fa {input[1]} results/{wildcards.sample}/binning/checkm_strainberry/
        """

rule sberry_bams:
    input:
        lathe_bam='results/{sample}/alignments/lathe-p1.bam',
        lathe_bamidx='results/{sample}/alignments/lathe-p1.bam.bai',
        bin_fasta='results/{sample}/binning/lathe-p1_bestbins/{bin}.fa',
        bin_faidx='results/{sample}/binning/lathe-p1_bestbins/{bin}.fa.fai',
    output:
        bin_bam='results/{sample}/alignments/lathe-p1_bestbins/{bin}.bam',
        bin_fastq='results/{sample}/fastqs/lathe-p1_bestbins/{bin}.fq.gz'
    conda:
        '../envs/alignment.yaml'
    threads:
        4
    shell:
        """
        mkdir -p results/{wildcards.sample}/alignments/lathe-p1_bestbins results/{wildcards.sample}/fastqs/lathe-p1_bestbins \
          && samtools view -bh -@{threads} -o {output.bin_bam} {input.lathe_bam} $(cut -f1 {input.bin_faidx} | tr "\\n" " ") \
          && samtools fastq {output.bin_bam} | gzip -c >{output.bin_fastq}
        """

rule sberry_separated_bins:
    input:
        fasta="results/{sample}/binning/lathe-p1_bestbins/{bin}.fa",
        faidx="results/{sample}/binning/lathe-p1_bestbins/{bin}.fa.fai",
        bam="results/{sample}/alignments/lathe-p1_bestbins/{bin}.bam",
        bamidx="results/{sample}/alignments/lathe-p1_bestbins/{bin}.bam.bai",
    output:
        sberry_fasta="results/{sample}/binning/strainberry/{bin}.fa",
        sberry_bam="results/{sample}/binning/strainberry/{bin}.bam",
    params:
        tech = '--nanopore' if config['technology'] == 'nanopore' else '',
    log:
        'logs/{sample}/strainberry/sberry_{bin}.log'
    benchmark:
        'benchmarks/{sample}/strainberry/sberry_{bin}.benchmark.txt'
    conda:
        "../envs/strainberry.yaml"
    threads:
        workflow.cores
    shell:
        """
        strainberry -r {input.fasta} -b {input.bam} -o results/{wildcards.sample}/binning/strainberry/sberry_{wildcards.bin} -c {threads} {params.tech} &>{log}
        if [ -f results/{wildcards.sample}/binning/strainberry/sberry_{wildcards.bin}/assembly.scaffolds.fa ]; then
          mv results/{wildcards.sample}/binning/strainberry/sberry_{wildcards.bin}/assembly.scaffolds.fa {output.sberry_fasta}
          mv results/{wildcards.sample}/binning/strainberry/sberry_{wildcards.bin}/assembly.scaffolds.fa.fai {output.sberry_fasta}.fai
          mv results/{wildcards.sample}/binning/strainberry/sberry_{wildcards.bin}/assembly.scaffolds.bam {output.sberry_bam}
          mv results/{wildcards.sample}/binning/strainberry/sberry_{wildcards.bin}/assembly.scaffolds.bam.bai {output.sberry_bam}.bai
        else
          cp {input.fasta} {output.sberry_fasta} && cp {input.faidx} {output.sberry_fasta}.fai
          cp {input.bam} {output.sberry_bam} && cp {input.bamidx} {output.sberry_bam}.bai
        fi
        """

rule sberry_bams_depth:
    input:
        sberry_bam='results/{sample}/binning/strainberry/{bin}.bam',
        sberry_fasta='results/{sample}/binning/strainberry/{bin}.fa',
        sberry_faidx='results/{sample}/binning/strainberry/{bin}.fa.fai',
    output:
        sberry_depth='results/{sample}/binning/strainberry/{bin}.depth',
    log:
        'logs/{sample}/strainberry/sberry_{bin}.depth.log'
    conda:
        "../envs/evaluation.yaml"
    threads:
        workflow.cores
    shell:
        """
        jgi_summarize_bam_contig_depths --minContigLength 1000 --minContigDepth 1 --percentIdentity 50 --outputDepth {output.sberry_depth} {input.sberry_bam} 2>{log}
        awk 'FNR==NR{{seq[$1]++;next}} seq[$1]>0{{print}}' {input.sberry_faidx} {output.sberry_depth} >{output.sberry_depth}.tmp
        mv {output.sberry_depth}.tmp {output.sberry_depth}
        """

rule sberry_classify_bins:
    input:
        fasta='results/{sample}/binning/strainberry/{bin}.fa',
    output:
        kraken2='results/{sample}/binning/strainberry/{bin}.kraken2',
    conda:
        "../envs/evaluation.yaml",
    threads:
        workflow.cores
    params:
        db=config['kraken2_db']
    shell:
        """
        kraken2 --db {params.db} --threads {threads} --output {output.kraken2} --use-names {input.fasta}
        """

def sberry_bin_depths(wildcards):
    checkpoint_output = checkpoints.sberry_bins_init.get(**wildcards).output.best_bins
    return expand('results/{sample}/binning/strainberry/{bin}.depth',
            sample=wildcards.sample,
            bin=glob_wildcards(os.path.join(checkpoint_output,'{binid}.fa')).binid)

def sberry_bin_class(wildcards):
    checkpoint_output = checkpoints.sberry_bins_init.get(**wildcards).output.best_bins
    return expand('results/{sample}/binning/strainberry/{bin}.kraken2',
            sample=wildcards.sample,
            bin=glob_wildcards(os.path.join(checkpoint_output,'{binid}.fa')).binid)

rule sberry_bins_stats:
    input:
        'results/{sample}/binning/strainberry',
        'results/{sample}/binning/strainberry.checkm.tsv',
        sberry_bin_depths,
        sberry_bin_class,
    output:
        bin_stats='results/{sample}/evaluation/strainberry.bin_stats.tsv',
    conda:
        "../envs/evaluation.yaml",
    shell:
        """
        mkdir -p results/{wildcards.sample}/evaluation \
          && python3 workflow/scripts/hsm_sberry_binstats.py --bins-dir {input[0]} --checkm {input[1]} --output {output.bin_stats}
        """

rule sberry_racon_init:
    input:  'results/{sample}/binning/strainberry/{bin}.fa',
    output: 'results/{sample}/binning/strainberry_racon/racon_0/{bin}.fa',
    shell:  'mkdir -p results/{wildcards.sample}/binning/strainberry_racon/racon_0 && ln -s "$(readlink -f {input})" {output}'

rule sberry_racon_paf:
    input:
        reads='results/{sample}/fastqs/lathe-p1_bestbins/{bin}.fq.gz',
        assembly='results/{sample}/binning/strainberry_racon/racon_{it}/{bin}.fa',
    output:
        pafgz='results/{sample}/binning/strainberry_racon/racon_{it}/{bin}.paf.gz'
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

def prev_racon_output(wildcards):
    it=int(wildcards.it)
    if it < 0:
        raise ValueError("racon iteration value should be >= 1")
    pafgz=f'results/{wildcards.sample}/binning/strainberry_racon/racon_{it-1}/{wildcards.bin}.paf.gz'
    assembly=f'results/{wildcards.sample}/binning/strainberry_racon/racon_{it-1}/{wildcards.bin}.fa'
    return (pafgz,assembly)

rule racon_consensus:
    input:
        'results/{sample}/fastqs/lathe-p1_bestbins/{bin}.fq.gz',
        prev_racon_output,
    output:
        'results/{sample}/binning/strainberry_racon/racon_{it}/{bin}.fa',
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

rule racon_alignment:
    input:
        fasta='results/{sample}/binning/strainberry_racon/racon_4/{bin}.fa',
        faidx='results/{sample}/binning/strainberry_racon/racon_4/{bin}.fa.fai',
        reads='results/{sample}/fastqs/lathe-p1_bestbins/{bin}.fq.gz',
    output:
        bam='results/{sample}/binning/strainberry_racon/racon_4/{bin}.bam',
    conda:
        '../envs/alignment.yaml'
    threads:
        workflow.cores
    params:
        preset = 'map-ont' if config['technology'] == 'nanopore' else 'map-pb',
    shell:
        """
        minimap2 -t {threads} -ax {params.preset} {input.fasta} {input.reads} | samtools sort --threads {threads} -o {output.bam}
        """

checkpoint medaka_init:
    input:
        fasta='results/{sample}/binning/strainberry_racon/racon_4/{bin}.fa',
        faidx='results/{sample}/binning/strainberry_racon/racon_4/{bin}.fa.fai',
    output:
        contigs=directory('results/{sample}/binning/strainberry_medaka/{bin}/contigs'),
        results=directory('results/{sample}/binning/strainberry_medaka/{bin}/results'),
        polished=directory('results/{sample}/binning/strainberry_medaka/{bin}/polished'),
    conda:
        '../envs/alignment.yaml'
    shell:
        """
        mkdir -p {output.contigs} {output.results} {output.polished} \
          && (cut -f1 {input.faidx} | xargs -n 1 -I{{}} bash -c 'samtools faidx {input.fasta} {{}} >{output.contigs}/{{}}.fa')
        """

rule medaka_consensus:
    input:
        'results/{sample}/binning/strainberry_medaka/{bin}/contigs/{contig}.fa',
        'results/{sample}/binning/strainberry_racon/racon_4/{bin}.bam',
        'results/{sample}/binning/strainberry_racon/racon_4/{bin}.bam.bai',
    output:
        'results/{sample}/binning/strainberry_medaka/{bin}/results/{contig}.hdf'
    conda:
        '../envs/polishing.yaml'
    shell:
        """
        medaka consensus {input[1]} '{output[0]}' --model r941_flip213 --threads 1 --regions '{wildcards.contig}:0-'
        """

rule medaka_stitch:
    input:  'results/{sample}/binning/strainberry_medaka/{bin}/results/{contig}.hdf',
    output: 'results/{sample}/binning/strainberry_medaka/{bin}/polished/{contig}.fa',
    conda:  '../envs/polishing.yaml'
    shell:  "medaka stitch '{input[0]}' '{output[0]}'"

def aggregate_medaka_results(wildcards):
    checkpoint_output = checkpoints.medaka_init.get(**wildcards).output.contigs
    return expand('results/{sample}/binning/strainberry_medaka/{bin}/polished/{contig}.fa',
            sample=wildcards.sample,
            bin=wildcards.bin,
            contig=glob_wildcards(os.path.join(checkpoint_output,'{ctg}.fa')).ctg)

rule medaka_aggregate:
    input:
        aggregate_medaka_results
    output:
        'results/{sample}/binning/strainberry_medaka/{bin}.fa'
    shell:
        """
        find results/{wildcards.sample}/binning/strainberry_medaka/{wildcards.bin}/polished -maxdepth 1 -type f -name '*.fa' -print0 \
          | xargs -0 cat | awk '/^>/{{ split($1,a,":"); printf("%s-%d\\n",a[1],++x[a[1]]); next }} {{print}}' >{output[0]}
        """

def sberry_medaka(wildcards):
    checkpoint_output = checkpoints.sberry_bins_init.get(**wildcards).output.best_bins
    return expand('results/{sample}/binning/strainberry_medaka/{bin}.fa',
            sample=wildcards.sample,
            bin=glob_wildcards(os.path.join(checkpoint_output,'{binid}.fa')).binid)

checkpoint sberry_medaka_done:
    input:  
        sberry_medaka
    output: 
        done='results/{sample}/binning/strainberry_medaka.done',
    shell:  
        'touch {output.done}'


rule sberry_medaka_alignment:
    input:
        fasta='results/{sample}/binning/strainberry_medaka/{bin}.fa',
        reads='results/{sample}/fastqs/lathe-p1_bestbins/{bin}.fq.gz',
    output:
        bam='results/{sample}/binning/strainberry_medaka/{bin}.bam',
    conda:
        '../envs/alignment.yaml'
    threads:
        workflow.cores
    params:
        preset = 'map-ont' if config['technology'] == 'nanopore' else 'map-pb',
    shell:
        """
        minimap2 -t {threads} -ax {params.preset} {input.fasta} {input.reads} | samtools sort --threads {threads} -o {output.bam}
        """

rule sberry_medaka_checkm_eval:
    input:
        'results/{sample}/binning/strainberry_medaka.done'
    output:
        checkm='results/{sample}/binning/strainberry_medaka.checkm.tsv'
    log:
        'logs/{sample}/strainberry_medaka_checkm.log'
    conda:
        "../envs/evaluation.yaml"
    threads:
        workflow.cores
    shell:
        """
        checkm lineage_wf --tab_table -f {output.checkm} -t {threads} --pplacer_threads {threads} \
          -x fa results/{wildcards.sample}/binning/strainberry_medaka/ results/{wildcards.sample}/binning/checkm_strainberry_medaka/
        """

rule sberry_medaka_depth:
    input:
        bam='results/{sample}/binning/strainberry_medaka/{bin}.bam',
        bamidx='results/{sample}/binning/strainberry_medaka/{bin}.bam.bai',
        fasta='results/{sample}/binning/strainberry_medaka/{bin}.fa',
        faidx='results/{sample}/binning/strainberry_medaka/{bin}.fa.fai',
    output:
        depth='results/{sample}/binning/strainberry_medaka/{bin}.depth',
    log:
        'logs/{sample}/strainberry_medaka/{bin}.depth.log'
    conda:
        "../envs/evaluation.yaml"
    threads:
        workflow.cores
    shell:
        """
        jgi_summarize_bam_contig_depths --minContigLength 1000 --minContigDepth 1 --percentIdentity 50 --outputDepth {output.depth} {input.bam} 2>{log}
        """

rule sberry_medaka_kraken:
    input:
        fasta='results/{sample}/binning/strainberry_medaka/{bin}.fa',
    output:
        kraken2='results/{sample}/binning/strainberry_medaka/{bin}.kraken2',
    conda:
        "../envs/evaluation.yaml",
    threads:
        workflow.cores
    params:
        db=config['kraken2_db']
    shell:
        """
        kraken2 --db {params.db} --threads {threads} --output {output.kraken2} --use-names {input.fasta}
        """

def sberry_medaka_depths(wildcards):
    checkpoint_output = checkpoints.sberry_medaka_done.get(**wildcards).output.done
    bin_ids=map( lambda path:os.path.splitext(os.path.basename(path))[0], glob.glob(f'results/{wildcards.sample}/binning/strainberry_medaka/*.fa') )
    return expand('results/{sample}/binning/strainberry_medaka/{bin}.depth', sample=wildcards.sample, bin=list(bin_ids))

def sberry_medaka_class(wildcards):
    checkpoint_output = checkpoints.sberry_medaka_done.get(**wildcards).output.done
    bin_ids=map( lambda path:os.path.splitext(os.path.basename(path))[0], glob.glob(f'results/{wildcards.sample}/binning/strainberry_medaka/*.fa') )
    return expand('results/{sample}/binning/strainberry_medaka/{bin}.kraken2', sample=wildcards.sample, bin=list(bin_ids))

rule sberry_medaka_binstats:
    input:
        'results/{sample}/binning/strainberry_medaka',
        'results/{sample}/binning/strainberry_medaka.checkm.tsv',
        sberry_medaka_depths,
        sberry_medaka_class,
    output:
        bin_stats='results/{sample}/evaluation/strainberry_medaka.bin_stats.tsv',
    conda:
        "../envs/evaluation.yaml",
    shell:
        """
        mkdir -p results/{wildcards.sample}/evaluation \
          && python3 workflow/scripts/hsm_sberry_binstats.py --bins-dir {input[0]} --checkm {input[1]} --output {output.bin_stats}
        """

rule sberry_racon_checkm_eval:
    input:
        'results/{sample}/binning/strainberry_medaka.done'
    output:
        checkm='results/{sample}/binning/strainberry_racon.checkm.tsv'
    log:
        'logs/{sample}/strainberry_racon_checkm.log'
    conda:
        "../envs/evaluation.yaml"
    threads:
        workflow.cores
    shell:
        """
        checkm lineage_wf --tab_table -f {output.checkm} -t {threads} --pplacer_threads {threads} \
          -x fa results/{wildcards.sample}/binning/strainberry_racon/racon_4/ results/{wildcards.sample}/binning/checkm_strainberry_racon/
        """

rule sberry_racon_depth:
    input:
        bam='results/{sample}/binning/strainberry_racon/racon_4/{bin}.bam',
        bamidx='results/{sample}/binning/strainberry_racon/racon_4/{bin}.bam.bai',
        fasta='results/{sample}/binning/strainberry_racon/racon_4/{bin}.fa',
        faidx='results/{sample}/binning/strainberry_racon/racon_4/{bin}.fa.fai',
    output:
        depth='results/{sample}/binning/strainberry_racon/racon_4/{bin}.depth',
    log:
        'logs/{sample}/strainberry_racon/{bin}.depth.log'
    conda:
        "../envs/evaluation.yaml"
    threads:
        workflow.cores
    shell:
        """
        jgi_summarize_bam_contig_depths --minContigLength 1000 --minContigDepth 1 --percentIdentity 50 --outputDepth {output.depth} {input.bam} 2>{log}
        """

rule sberry_racon_kraken:
    input:
        fasta='results/{sample}/binning/strainberry_racon/racon_4/{bin}.fa',
    output:
        kraken2='results/{sample}/binning/strainberry_racon/racon_4/{bin}.kraken2',
    conda:
        "../envs/evaluation.yaml",
    threads:
        workflow.cores
    params:
        db=config['kraken2_db']
    shell:
        """
        kraken2 --db {params.db} --threads {threads} --output {output.kraken2} --use-names {input.fasta}
        """

def sberry_racon_depths(wildcards):
    checkpoint_output = checkpoints.sberry_medaka_done.get(**wildcards).output.done
    bin_ids=map( lambda path:os.path.splitext(os.path.basename(path))[0], glob.glob(f'results/{wildcards.sample}/binning/strainberry_racon/racon_4/*.fa') )
    return expand('results/{sample}/binning/strainberry_racon/racon_4/{bin}.depth', sample=wildcards.sample, bin=list(bin_ids))

def sberry_racon_class(wildcards):
    checkpoint_output = checkpoints.sberry_medaka_done.get(**wildcards).output.done
    bin_ids=map( lambda path:os.path.splitext(os.path.basename(path))[0], glob.glob(f'results/{wildcards.sample}/binning/strainberry_racon/racon_4/*.fa') )
    return expand('results/{sample}/binning/strainberry_racon/racon_4/{bin}.kraken2', sample=wildcards.sample, bin=list(bin_ids))

rule sberry_racon_binstats:
    input:
        'results/{sample}/binning/strainberry_racon/racon_4',
        'results/{sample}/binning/strainberry_racon.checkm.tsv',
        sberry_racon_depths,
        sberry_racon_class,
    output:
        bin_stats='results/{sample}/evaluation/strainberry_racon.bin_stats.tsv',
    conda:
        "../envs/evaluation.yaml",
    shell:
        """
        mkdir -p results/{wildcards.sample}/evaluation \
          && python3 workflow/scripts/hsm_sberry_binstats.py --bins-dir {input[0]} --checkm {input[1]} --output {output.bin_stats}
        """



