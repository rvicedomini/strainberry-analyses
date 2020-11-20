rule mock3_barplots:
    input:
        "results/mock3/assembly_eval/canu.report.tsv",
        "results/mock3/assembly_eval/sberry_canu_n2_scf.report.tsv",
        "results/mock3/assembly_eval/flye.report.tsv",
        "results/mock3/assembly_eval/sberry_flye_n2_scf.report.tsv",
    output:
        "results/mock3/assembly_eval/mock3.asmstats.tsv",
        "results/mock3/assembly_eval/mock3.ani.barplot.pdf",
        "results/mock3/assembly_eval/mock3.dupratio.barplot.pdf",
    conda:
        "../envs/assembly_eval.yaml"
    shell:
        """
        for f in {input}; do
            fbase=${{f##*/}}
            asmname=${{fbase%%.*}}
            egrep -v '^#' ${{f}} | egrep 'Ecoli-K12|Ecoli-W' | awk -v asmname=${{asmname}} '{{OFS="\\t";print asmname,$0}}'
        done >{output[0]}
        python3 workflow/scripts/assembly_barplots.py -i {output[0]} -p results/mock3/assembly_eval/mock3
        """

rule mock3_circos_config:
    input:  
        'resources/mock/circos/etc/mock3.conf',
        'resources/mock/circos/data/mock3.karyotype.txt',
    output: 
        'results/mock3/assembly_eval/circos/etc/mock3.conf',
        'results/mock3/assembly_eval/circos/data/mock3.karyotype.txt'
    shell:  
        """
        mkdir -p results/mock3/assembly_eval/circos \
          && cp -rn resources/mock/circos/{{data,etc}} results/mock3/assembly_eval/circos
        """

rule mock3_circos_prepare:
    input:
        'results/mock3/assembly_eval/circos/etc/mock3.conf',
        'results/mock3/assembly_eval/{assembly}.report.tsv',
    output:
        "results/mock3/assembly_eval/circos/data/{assembly}.1coords.tsv",
        "results/mock3/assembly_eval/circos/data/{assembly}.snps.int.tsv",
    shell:
        """
        : > {output[0]}
        find results/mock3/assembly_eval/{wildcards.assembly}/dnadiff/ -maxdepth 1 -name '*.1coords' -exec cat {{}} + \
          | awk 'BEGIN{{OFS="\\t"}} $2!="." && $3!="."{{ print $12,$1,$2,$7 }}' \
          | tee --append {output[0]} >/dev/null
        : > {output[1]}
        for f in $(find results/mock3/assembly_eval/{wildcards.assembly}/dnadiff/ -maxdepth 1 -name '*.snps'); do
          cat "${{f}}" | sort -k11,11 -k1,1n \
            | awk 'BEGIN{{beg=1}} $2!="." && $3!="."{{ step=int($7/1500)+1;end=beg+step;while($1>end){{beg=end+1;end=beg+step}} snps[beg]+=1;x[beg]=end}} END{{ for(beg in snps){{printf("%s\\t%d\\t%d\\t%d\\n",$11,beg,x[beg],snps[beg])}} }}'
        done | sort -k1,1 -k2,2n | tee --append {output[1]} >/dev/null
        """

rule mock3_circos_plot:
    input: 
        'results/mock3/assembly_eval/circos/etc/mock3.conf',
        # aligned blocks
        'results/mock3/assembly_eval/circos/data/canu.1coords.tsv',
        'results/mock3/assembly_eval/circos/data/sberry_canu_n2_scf.1coords.tsv',
        'results/mock3/assembly_eval/circos/data/flye.1coords.tsv',
        'results/mock3/assembly_eval/circos/data/sberry_flye_n2_scf.1coords.tsv',
        # snps intervals
        'results/mock3/assembly_eval/circos/data/canu.snps.int.tsv',
        'results/mock3/assembly_eval/circos/data/sberry_canu_n2_scf.snps.int.tsv',
        'results/mock3/assembly_eval/circos/data/flye.snps.int.tsv',
        'results/mock3/assembly_eval/circos/data/sberry_flye_n2_scf.snps.int.tsv',
    output: 
        'results/mock3/assembly_eval/mock3_circos.svg'
    log:
        'logs/mock9/circos_mock3.log'
    conda: 
        '../envs/circos.yaml'
    shell: 
        """
        cd results/mock3/assembly_eval/circos \
          && env PERL5LIB= PERL_LOCAL_LIB_ROOT= circos -conf etc/mock3.conf >{log} \
          && cp mock3_circos.svg ../
        """

