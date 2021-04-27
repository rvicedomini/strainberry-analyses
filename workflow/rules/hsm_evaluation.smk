
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
        lathe='results/{sample}/evaluation/lathe-p1.bin_stats.tsv',
        sberry='results/{sample}/evaluation/strainberry_medaka.bin_stats.tsv',
    output:
        'results/{sample}/evaluation/lathe-p1_strainberry.barplot.pdf',
        'results/{sample}/evaluation/lathe-p1_strainberry.barplot.svg',
    conda:
        "../envs/evaluation.yaml"
    shell:
        """
        python3 workflow/scripts/hsm_barplot.py --pre {input.lathe} --sep {input.sberry} --prefix results/{wildcards.sample}/evaluation/lathe-p1_strainberry.barplot
        """

rule hsm_vatypica:
    input:
        lathe_binstats='results/{sample}/evaluation/lathe-p1.bin_stats.tsv',
        lathe_seq2bin='results/{sample}/binning/lathe-p1.bins.tsv',
        lathe_depth='results/{sample}/binning/lathe-p1.depth.txt',
        sberry_done='results/{sample}/binning/strainberry_medaka.done',
    output:
        outdir=directory('results/{sample}/vatypica'),
        plot='results/{sample}/evaluation/strainberry_lathe-p1.vatypica.depth.pdf'
    conda:
        "../envs/evaluation.yaml"
    shell:
        """
        mkdir -p {output.outdir}
        vatypica_bin="$(grep -m1 -F 'Veillonella_atypica' {input.lathe_binstats} | awk '{{printf("%s",$1);exit}}')"

        awk -v binid="${{vatypica_bin}}" '$2==binid{{print $1}}' {input.lathe_seq2bin} > {output.outdir}/lathe-p1.vatypica.list
        awk 'BEGIN{{OFS="\\t"}} FNR==NR{{seq[$1]++;next}} ($1 in seq){{print $1,$2,$3}}' {output.outdir}/lathe-p1.vatypica.list {input.lathe_depth} >{output.outdir}/lathe-p1.vatypica.depth.tsv
        
        cut -f1 results/{wildcards.sample}/binning/strainberry_medaka/"${{vatypica_bin}}".fa.fai > {output.outdir}/strainberry_medaka.vatypica.list
        cut -f1,2,3 results/{wildcards.sample}/binning/strainberry_medaka/"${{vatypica_bin}}".depth | tail -q -n+2 >{output.outdir}/strainberry_medaka.vatypica.depth.tsv

        python3 workflow/scripts/hsm_plot_coverage.py --before {output.outdir}/lathe-p1.vatypica.depth.tsv --after {output.outdir}/strainberry_medaka.vatypica.depth.tsv --prefix results/{wildcards.sample}/evaluation/strainberry_lathe-p1.vatypica.depth
        """

rule hsm_vatypica_mummerplot:
    input:
        lathe_binstats='results/{sample}/evaluation/lathe-p1.bin_stats.tsv',
        strainberry_binstats='results/{sample}/evaluation/strainberry_medaka.bin_stats.tsv',
        vasch6_ref='resources/{sample}/references/NZ_AEDR00000000.1.fasta',
        vacol7_ref='resources/{sample}/references/NZ_AEDS00000000.1.fasta'
    output:
        outdir=directory('results/{sample}/vatypica'),
        vasch6='results/{sample}/evaluation/strainberry_medaka.vatypica.ACS-049-V-Sch6.ps',
        vacol7='results/{sample}/evaluation/strainberry_medaka.vatypica.ACS-134-V-Col7a.ps'
    conda:
        "../envs/mummerplot.yaml"
    shell:
        """
        mkdir -p {output.outdir}
        vatypica_bin="$(grep -m1 -F 'Veillonella_atypica' {input.lathe_binstats} | awk '{{printf("%s",$1);exit}}')"

        grep -F 'ACS-049-V-Sch6' results/{wildcards.sample}/binning/strainberry_medaka/"${{vatypica_bin}}".kraken2 | cut -f2 | xargs samtools faidx results/{wildcards.sample}/binning/strainberry_medaka/"${{vatypica_bin}}".fa > {output.outdir}/strainberry_medaka.ACS-049-V-Sch6.fa
        nucmer --maxmatch -p {output.outdir}/strainberry_medaka.ACS-049-V-Sch6 {input.vasch6_ref} {output.outdir}/strainberry_medaka.ACS-049-V-Sch6.fa
        mummerplot -c --layout --postscript -p {output.outdir}/strainberry_medaka.ACS-049-V-Sch6 -R {input.vasch6_ref} -Q {output.outdir}/strainberry_medaka.ACS-049-V-Sch6.fa {output.outdir}/strainberry_medaka.ACS-049-V-Sch6.delta
        mv {output.outdir}/strainberry_medaka.ACS-049-V-Sch6.ps {output.vasch6}

        grep -F 'ACS-134-V-Col7a' results/{wildcards.sample}/binning/strainberry_medaka/"${{vatypica_bin}}".kraken2 | cut -f2 | xargs samtools faidx results/{wildcards.sample}/binning/strainberry_medaka/"${{vatypica_bin}}".fa > {output.outdir}/strainberry_medaka.ACS-134-V-Col7a.fa
        nucmer --maxmatch -p {output.outdir}/strainberry_medaka.ACS-134-V-Col7a {input.vacol7_ref} {output.outdir}/strainberry_medaka.ACS-134-V-Col7a.fa
        mummerplot -c --layout --postscript -p {output.outdir}/strainberry_medaka.ACS-134-V-Col7a -R {input.vacol7_ref} -Q {output.outdir}/strainberry_medaka.ACS-134-V-Col7a.fa {output.outdir}/strainberry_medaka.ACS-134-V-Col7a.delta
        mv {output.outdir}/strainberry_medaka.ACS-134-V-Col7a.ps {output.vacol7}
        """


rule hsm_eeligens:
    input:
        lathe_binstats='results/{sample}/evaluation/lathe-p1.bin_stats.tsv',
        lathe_seq2bin='results/{sample}/binning/lathe-p1.bins.tsv',
        lathe_depth='results/{sample}/binning/lathe-p1.depth.txt',
        sberry_done='results/{sample}/binning/strainberry_medaka.done',
    output:
        outdir=directory('results/{sample}/eeligens'),
        plot='results/{sample}/evaluation/strainberry_lathe-p1.eeligens.depth.pdf'
    conda:
        "../envs/evaluation.yaml"
    shell:
        """
        mkdir -p {output.outdir}
        eeligens_bin="$(grep -m1 -F '[Eubacterium]_eligens' {input.lathe_binstats} | awk '{{printf("%s",$1);exit}}')"

        awk -v binid="${{eeligens_bin}}" '$2==binid{{print $1}}' {input.lathe_seq2bin} > {output.outdir}/lathe-p1.eeligens.list
        awk 'BEGIN{{OFS="\\t"}} FNR==NR{{seq[$1]++;next}} ($1 in seq){{print $1,$2,$3}}' {output.outdir}/lathe-p1.eeligens.list {input.lathe_depth} >{output.outdir}/lathe-p1.eeligens.depth.tsv

        cut -f1 results/{wildcards.sample}/binning/strainberry_medaka/"${{eeligens_bin}}".fa.fai > {output.outdir}/strainberry_medaka.eeligens.list
        cut -f1,2,3 results/{wildcards.sample}/binning/strainberry_medaka/"${{eeligens_bin}}".depth | tail -q -n+2 >{output.outdir}/strainberry_medaka.eeligens.depth.tsv
        
        python3 workflow/scripts/hsm_plot_coverage.py --before {output.outdir}/lathe-p1.eeligens.depth.tsv --after {output.outdir}/strainberry_medaka.eeligens.depth.tsv --prefix results/{wildcards.sample}/evaluation/strainberry_lathe-p1.eeligens.depth
        """


