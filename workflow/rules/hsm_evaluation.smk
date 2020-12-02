
rule hsm_kraken2:
    input:
        fasta='results/{sample}/assemblies/{assembly}.fa',
    output:
        kreport='results/{sample}/kraken2/{assembly}.kreport',
        kraken2='results/{sample}/kraken2/{assembly}.kraken2',
    conda:
        '../envs/evaluation.yaml',
    threads:
        workflow.cores
    params:
        db=config['kraken2_db']
    shell:
        """
        kraken2 --db {params.db} --threads {threads} --report {output.kreport} --output {output.kraken2} --use-names {input.fasta}
        """

rule hsm_lathe_binning:
    input:
        fasta='results/hsm/assemblies/lathe-p1.fa',
        faidx='results/hsm/assemblies/lathe-p1.fa.fai',
        bam='results/hsm/alignments/lathe-p1.bam',
        bamidx='results/hsm/alignments/lathe-p1.bam.bai',
    output:
        depth='results/hsm/binning/lathe-p1.depth.txt',
        bins=directory('results/hsm/binning/lathe-p1'),
        seq2bin='results/hsm/binning/lathe-p1.bins.tsv',
    log:
        'logs/hsm/lathe-p1_binning.log'
    conda:
        "../envs/evaluation.yaml"
    threads:
        workflow.cores
    shell:
        """
        mkdir -p {output.bins} \
          && jgi_summarize_bam_contig_depths --minContigLength 1000 --minContigDepth 1 --percentIdentity 50 --outputDepth {output.depth} {input.bam} 2>{log} \
          && metabat2 -t {threads} -i {input.fasta} -a {output.depth} -o {output.bins}/bin --unbinned --saveCls --seed 1 2>>{log} \
          && python3 workflow/scripts/seq2bin.py -i {output.bins} -x fa -o {output.seq2bin} 2>>{log}
        """

# Strainberry bins are not created with metabat2
# In order to compare species/strain assemblies before/after the separation,
# Strainberry bins generated with the separated sequences of each Lathe bin
rule hsm_sberry_binning:
    input:
        sberry_fasta='results/hsm/assemblies/sberry_lathe-p1_n2_ctg.medaka.fa',
        bam='results/hsm/alignments/sberry_lathe-p1_n2_ctg.medaka.bam',
        bamidx='results/hsm/alignments/sberry_lathe-p1_n2_ctg.medaka.bam.bai',
        lathe_seq2bin='results/hsm/binning/lathe-p1.bins.tsv',
    output:
        depth='results/hsm/binning/sberry_lathe-p1_n2_ctg.medaka.depth.txt',
        bins=directory('results/hsm/binning/sberry_lathe-p1_n2_ctg.medaka'),
        seq2bin='results/hsm/binning/sberry_lathe-p1_n2_ctg.medaka.bins.tsv',
    log:
        'logs/hsm/sberry_lathe-p1_n2_ctg.medaka_binning.log'
    conda:
        "../envs/evaluation.yaml"
    threads:
        workflow.cores
    shell:
        """
        jgi_summarize_bam_contig_depths --minContigLength 1000 --minContigDepth 1 --percentIdentity 50 --outputDepth {output.depth} {input.bam} 2>{log} \
          && python3 workflow/scripts/sberry_map2bins.py -i {input.sberry_fasta} -s {input.lathe_seq2bin} -b {output.bins} -x fa -f 2>>{log} \
          && python3 workflow/scripts/seq2bin.py -i {output.bins} -x fa -o {output.seq2bin} 2>>{log}
        """

rule hsm_checkm_eval:
    input:
        depth='results/{sample}/binning/{assembly}.depth.txt',
        bins='results/{sample}/binning/{assembly}',
    output:
        checkm='results/{sample}/binning/{assembly}.checkm.tsv'
    log:
        'logs/{sample}/{assembly}_checkm.log'
    conda:
        "../envs/evaluation.yaml"
    threads:
        workflow.cores
    shell:
        """
        checkm lineage_wf --tab_table -f {output.checkm} -t {threads} --pplacer_threads {threads} \
          -x fa {input.bins} results/{wildcards.sample}/binning/checkm_{wildcards.assembly}/
        """

rule hsm_classify_bins:
    input:
        kraken2='results/{sample}/kraken2/{assembly}.kraken2',
        seq2bin='results/{sample}/binning/{assembly}.bins.tsv',
    output:
        classif='results/{sample}/binning/{assembly}.class.tsv'
    conda:
        "../envs/evaluation.yaml"
    shell:
        """
        python3 workflow/scripts/classify_bin.py --kraken {input.kraken2} --seq2bin {input.seq2bin} -o {output.classif} \
          && LC_ALL=C sort -k1,1V -k7,7gr -o {output.classif} {output.classif}
        """

rule hsm_bin_stats:
    input:
        depth='results/{sample}/binning/{assembly}.depth.txt',
        seq2bin='results/{sample}/binning/{assembly}.bins.tsv',
        checkm='results/{sample}/binning/{assembly}.checkm.tsv',
        classif='results/{sample}/binning/{assembly}.class.tsv',
    output:
        bin_stats='results/{sample}/evaluation/{assembly}.bin_stats.tsv',
    params:
        keep_best_only = lambda w: '--keep-best-only' if w.assembly=='lathe-p1' else ''
    conda:
        "../envs/evaluation.yaml"
    shell:
        """
        mkdir -p results/{wildcards.sample}/evaluation \
          && python3 workflow/scripts/bin_stats.py {params.keep_best_only} \
               --depth {input.depth} --seq2bin {input.seq2bin} --checkm {input.checkm} --class {input.classif} \
               --output {output.bin_stats}
        """

rule hsm_barplot:
    input:
        lathe='results/hsm/evaluation/lathe-p1.bin_stats.tsv',
        sberry='results/hsm/evaluation/sberry_lathe-p1_n2_ctg.medaka.bin_stats.tsv',
    output:
        'results/hsm/evaluation/sberry_lathe-p1.barplot.pdf',
        'results/hsm/evaluation/sberry_lathe-p1.barplot.svg',
    conda:
        "../envs/evaluation.yaml"
    shell:
        """
        python3 workflow/scripts/hsm_barplot.py --pre {input.lathe} --sep {input.sberry} --prefix results/hsm/evaluation/sberry_lathe-p1.barplot
        """

rule hsm_vatypica:
    input:
        lathe_binstats='results/hsm/evaluation/lathe-p1.bin_stats.tsv',
        lathe_seq2bin='results/hsm/binning/lathe-p1.bins.tsv',
        lathe_depth='results/hsm/binning/lathe-p1.depth.txt',
        sberry_seq2bin='results/hsm/binning/sberry_lathe-p1_n2_ctg.medaka.bins.tsv',
        sberry_depth='results/hsm/binning/sberry_lathe-p1_n2_ctg.medaka.depth.txt',
    output:
        outdir=directory('results/hsm/vatypica'),
        plot='results/hsm/evaluation/sberry_lathe-p1.vatypica.depth.pdf'
    conda:
        "../envs/evaluation.yaml"
    shell:
        """
        mkdir -p {output.outdir}
        vatypica_bin="$(fgrep -m1 'Veillonella_atypica' {input.lathe_binstats} | awk '{{print $1;exit}}')"

        awk -v binid="${{vatypica_bin}}" '$2==binid{{print $1}}' {input.lathe_seq2bin} > {output.outdir}/lathe-p1.vatypica.list
        awk 'BEGIN{{OFS="\\t"}} FNR==NR{{seq[$1]++;next}} ($1 in seq){{print $1,$2,$3}}' {output.outdir}/lathe-p1.vatypica.list {input.lathe_depth} >{output.outdir}/lathe-p1.vatypica.depth.tsv

        awk -v binid="${{vatypica_bin}}" '$2==binid{{print $1}}' {input.sberry_seq2bin} > {output.outdir}/sberry_lathe-p1.vatypica.list
        awk 'BEGIN{{OFS="\\t"}} FNR==NR{{seq[$1]++;next}} ($1 in seq){{print $1,$2,$3}}' {output.outdir}/sberry_lathe-p1.vatypica.list {input.sberry_depth} >{output.outdir}/sberry_lathe-p1.vatypica.depth.tsv

        python3 workflow/scripts/plot_coverage_diff.py --before {output.outdir}/lathe-p1.vatypica.depth.tsv --after {output.outdir}/sberry_lathe-p1.vatypica.depth.tsv --prefix results/hsm/evaluation/sberry_lathe-p1.vatypica.depth
        """

rule hsm_eeligens:
    input:
        lathe_binstats='results/hsm/evaluation/lathe-p1.bin_stats.tsv',
        lathe_seq2bin='results/hsm/binning/lathe-p1.bins.tsv',
        lathe_depth='results/hsm/binning/lathe-p1.depth.txt',
        sberry_seq2bin='results/hsm/binning/sberry_lathe-p1_n2_ctg.medaka.bins.tsv',
        sberry_depth='results/hsm/binning/sberry_lathe-p1_n2_ctg.medaka.depth.txt',
    output:
        outdir=directory('results/hsm/eeligens'),
        plot='results/hsm/evaluation/sberry_lathe-p1.eeligens.depth.pdf'
    conda:
        "../envs/evaluation.yaml"
    shell:
        """
        mkdir -p {output.outdir}
        eeligens_bin="$(fgrep -m1 '[Eubacterium]_eligens' {input.lathe_binstats} | awk '{{print $1;exit}}')"

        awk -v binid="${{eeligens_bin}}" '$2==binid{{print $1}}' {input.lathe_seq2bin} > {output.outdir}/lathe-p1.eeligens.list
        awk 'BEGIN{{OFS="\\t"}} FNR==NR{{seq[$1]++;next}} ($1 in seq){{print $1,$2,$3}}' {output.outdir}/lathe-p1.eeligens.list {input.lathe_depth} >{output.outdir}/lathe-p1.eeligens.depth.tsv

        awk -v binid="${{eeligens_bin}}" '$2==binid{{print $1}}' {input.sberry_seq2bin} > {output.outdir}/sberry_lathe-p1.eeligens.list
        awk 'BEGIN{{OFS="\\t"}} FNR==NR{{seq[$1]++;next}} ($1 in seq){{print $1,$2,$3}}' {output.outdir}/sberry_lathe-p1.eeligens.list {input.sberry_depth} >{output.outdir}/sberry_lathe-p1.eeligens.depth.tsv

        python3 workflow/scripts/plot_coverage_diff.py --before {output.outdir}/lathe-p1.eeligens.depth.tsv --after {output.outdir}/sberry_lathe-p1.eeligens.depth.tsv --prefix results/hsm/evaluation/sberry_lathe-p1.eeligens.depth
        """


