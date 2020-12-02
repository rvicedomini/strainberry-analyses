rule mock9_n2_barplots:
    input:
        "results/mock9/evaluation/canu.report.tsv",
        "results/mock9/evaluation/sberry_canu_n2_scf.report.tsv",
        "results/mock9/evaluation/flye.report.tsv",
        "results/mock9/evaluation/sberry_flye_n2_scf.report.tsv",
    output:
        "results/mock9/evaluation/mock9_n2.asmstats.tsv",
        "results/mock9/evaluation/mock9_n2.ani.barplot.pdf",
        "results/mock9/evaluation/mock9_n2.dupratio.barplot.pdf",
    conda:
        "../envs/evaluation.yaml"
    shell:
        """
        for f in {input}; do
            fbase=${{f##*/}}
            asmname=${{fbase%%.*}}
            egrep -v '^#' ${{f}} | egrep 'Kpneumoniae|Saureus-ATCC|Saureus-FDAA' \
              | awk -v asmname=${{asmname}} '{{OFS="\\t";print asmname,$0}}'
        done >{output[0]}
        python3 workflow/scripts/assembly_barplots.py -i {output[0]} -p results/mock9/evaluation/mock9_n2
        """

rule mock9_n2_circos_config:
    input:  
        'resources/mock/circos/etc/mock9_n2.conf',
        'resources/mock/circos/data/mock9_n2.karyotype.txt',
    output: 
        'results/mock9/evaluation/circos_mock9_n2/etc/mock9_n2.conf',
        'results/mock9/evaluation/circos_mock9_n2/data/mock9_n2.karyotype.txt',
    shell:
        """
        mkdir -p results/mock9/evaluation/circos_mock9_n2 \
          && cp -rn resources/mock/circos/{{data,etc}} results/mock9/evaluation/circos_mock9_n2
        """

rule mock9_n2_circos_prepare:
    input:
        'results/mock9/evaluation/{assembly}.report.tsv',
    output:
        "results/mock9/evaluation/circos_mock9_n2/data/{assembly}.1coords.tsv",
        "results/mock9/evaluation/circos_mock9_n2/data/{assembly}.snps.int.tsv",
    shell:
        """
        : > {output[0]}
        find results/mock9/evaluation/{wildcards.assembly}/dnadiff/ -maxdepth 1 -name '*.1coords' -exec cat {{}} + \
          | awk 'BEGIN{{OFS="\\t"}} $2!="." && $3!="."{{ print $12,$1,$2,$7 }}' \
          | tee --append {output[0]} >/dev/null
        : > {output[1]}
        for f in $(find results/mock9/evaluation/{wildcards.assembly}/dnadiff/ -maxdepth 1 -name '*.snps'); do
          cat "${{f}}" | sort -k11,11 -k1,1n \
            | awk 'BEGIN{{beg=1}} $2!="." && $3!="."{{ step=int($7/1500)+1;end=beg+step;while($1>end){{beg=end+1;end=beg+step}} snps[beg]+=1;x[beg]=end}} END{{ for(beg in snps){{printf("%s\\t%d\\t%d\\t%d\\n",$11,beg,x[beg],snps[beg])}} }}'
        done | sort -k1,1 -k2,2n | tee --append {output[1]} >/dev/null
        """

rule mock9_n2_circos_plot:
    input:
        'results/mock9/evaluation/circos_mock9_n2/etc/mock9_n2.conf',
        # aligned blocks
        'results/mock9/evaluation/circos_mock9_n2/data/canu.1coords.tsv',
        'results/mock9/evaluation/circos_mock9_n2/data/sberry_canu_n2_scf.1coords.tsv',
        'results/mock9/evaluation/circos_mock9_n2/data/flye.1coords.tsv',
        'results/mock9/evaluation/circos_mock9_n2/data/sberry_flye_n2_scf.1coords.tsv',
        # snps intervals
        'results/mock9/evaluation/circos_mock9_n2/data/canu.snps.int.tsv',
        'results/mock9/evaluation/circos_mock9_n2/data/sberry_canu_n2_scf.snps.int.tsv',
        'results/mock9/evaluation/circos_mock9_n2/data/flye.snps.int.tsv',
        'results/mock9/evaluation/circos_mock9_n2/data/sberry_flye_n2_scf.snps.int.tsv',
    output: 
        'results/mock9/evaluation/mock9_n2_circos.svg'
    log:
        'logs/mock9/circos_mock9_n2.log'
    conda: 
        '../envs/circos.yaml'
    shell: 
        """
        logfile="$(readlink -f {log})"
        cd results/mock9/evaluation/circos_mock9_n2 \
          && env PERL5LIB= PERL_LOCAL_LIB_ROOT= circos -conf etc/mock9_n2.conf &>"${{logfile}}" \
          && cp mock9_n2_circos.svg ../
        """




