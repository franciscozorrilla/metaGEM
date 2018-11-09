configfile: "config.yaml"

IDs = "Sample118_s1e5_R1.fa Sample120_s1e5_R1.fa Sample127_s1e5_R1.fa Sample134_s1e5_R1.fa Sample177_s1e5_R1.fa Sample215_s1e5_R1.fa Sample230_s1e5_R1.fa Sample234_s1e5_R1.fa Sample244_s1e5_R1.fa Sample261_s1e5_R1.fa Sample263_s1e5_R1.fa Sample290_s1e5_R1.fa Sample302_s1e5_R1.fa Sample321_s1e5_R1.fa Sample330_s1e5_R1.fa Sample343_s1e5_R1.fa".split()

#import os
#import glob
#SPECIES = glob.glob('carvemeOut/*.txt')
#SPECIES = [os.path.splitext(val)[0] for val in glob.glob('carvemeOut/*.txt')]

rule all:
    input: "fetchMG/output", dynamic("carvemeOut/{species}.xml.html")
    shell:"snakemake --dag | dot -Tpng > pipemap.png"

rule mergeReads:
    output:"All_{read}.fa"
    shell:"cat {config[paths][raw_reads]}/Sample*_{wildcards.read}.fa  > {output}"

rule megahit:
    input:R1="All_R1.fa",R2="All_R2.fa"
    output:"megahitAssembly"
    shell:"megahit --min-count {config[megahit_params][mincount]} --k-list {config[megahit_params][klistdefault]} -1 {input.R1} -2 {input.R2} -o {output}"

rule cutcontigs:
    input:"megahitAssembly"
    output:"megahit_c10K.fa"
    shell:"set +u;source activate concoct_env;set -u; cd {input}; python {config[cutcontigs_params][script_dir]} -c {config[cutcontigs_params][chunk_size]} -o {config[cutcontigs_params][o]} -m final.contigs.fa > {output}; cp {output} {config[paths][concoct_run]}; cd .."

rule bowtie:
    input:"megahit_c10K.fa"
    output:"map/{id}"
    shell:
        """
        export MRKDUP={config[bowtie_params][MRKDUP_jardir]};
        set +u;source activate concoct_env;set -u;
        bowtie2-build {input} {input};
        mkdir -p {output}
        cd {output};
        bash {config[bowtie_params][MRKDUP_shelldir]} -ct 1 -p '-f' {config[paths][raw_reads]}/$(basename /{output}) $(echo {config[paths][raw_reads]}/$(basename /{output}) | sed s/R1/R2/) pair {config[paths][concoct_run]}/{input} asm bowtie2;
        cd {config[paths][concoct_run]};
        """

rule covtable:
    input: expand("map/{id}", id=IDs)
    output:"concoct-input/concoct_inputtable.tsv"
    shell:
        """
        set +u;source activate concoct_env;set -u;
        cd {config[paths][concoct_run]}/{config[bowtie_params][outputdir]}
        python {config[paths][CONCOCT]}/scripts/gen_input_table.py --isbedfiles \
            --samplenames <(for s in Sample*; do echo $s | cut -d'_' -f1; done) \
            ../megahit_c10K.fa */bowtie2/asm_pair-smds.coverage > concoct_inputtable.tsv;
        mkdir -p {config[paths][concoct_run]}/{config[covtable_params][outdir]};
        mv concoct_inputtable.tsv {config[paths][concoct_run]}/{config[covtable_params][outdir]}/;
        """

#rule linktable:
#    input:"contigs/velvet_71_c10K.fa"
#    output:"concoct-input/concoct_linkage.tsv"
#    shell:
#        """
#        set +u;source activate concoct_env;set -u;
#        cd {config[paths][concoct_run]}/map;
#        python {config[paths][CONCOCT]}/scripts/bam_to_linkage.py -m {config[linktable_params][m]} \
#        --regionlength {config[linktable_params][regionlength]} --{config[linktable_params][search]} \
#        --samplenames <(for s in Sample*; do echo $s | cut -d'_' -f1; done) \
#        ../{input} Sample*/bowtie2/asm_pair-smds.bam > concoct_linkage.tsv;
#        mv concoct_linkage.tsv {config[paths][concoct_run]}/concoct-input/;
#        """

rule concoct:
    input:"concoct-input/concoct_inputtable.tsv"
    output:"concoct-output/clustering_gt1000.csv"
    shell:
        """
        set +u;source activate concoct_env;set -u;
        cut {config[concoct_params][cutoptions]} {input} > {config[concoct_params][covfile]};
        concoct -c {config[concoct_params][c]} --coverage_file {config[concoct_params][covfile]} --composition_file {config[concoct_params][compositionfile]} -b {config[concoct_params][outputdir]}
        """

rule evaluate:
    input:"concoct-output/clustering_gt1000.csv"
    output:"evaluation-output/ClusterPlot.pdf"
    shell:
        """
        set +u;source activate concoct_env;set -u
        mkdir -p {config[evaluate_params][outputdir]};
        Rscript {config[paths][CONCOCT]}/scripts/ClusterPlot.R -c {input} -p {config[evaluate_params][p]} -m {config[evaluate_params][m]} -r {config[evaluate_params][r]} -l -o {output};
        cp {config[paths][CONCOCT_TEST]}/evaluation-output/clustering_gt1000_s.csv {config[evaluate_params][outputdir]}/;
        """
        #NOTE: this line throws cryptic error {config[paths][CONCOCT]}/scripts/Validate.pl --cfile={config[evaluate_params][cfile]} --sfile={config[evaluate_params][sfile]} --ofile={config[evaluate_params][ofile]} --ffile={config[evaluate_params][ffile]};
        #NOTE: this evaluation is based on a toy example where the species are already know, need to replace this with TAXAssign or similar tool for data where species are not known
        #NOTE: heatmap doesnt work:  Rscript {config[paths][CONCOCT]}/scripts/ConfPlot.R  -c evaluation-output/clustering_gt1000_conf.csv -o  evaluation-output/clustering_gt1000_conf.pdf;

rule prodigal:
    input:"megahit_c10K.fa"
    output:"annotations/proteins/megahit_c10K.gff"
    shell:
        """
        set +u;source activate concoct_env;set -u
        cd {config[paths][concoct_run]}
        mkdir -p {config[paths][concoct_run]}/{config[prodigal_params][outputdir]}
        prodigal -a {config[prodigal_params][outputdir]}/megahit_c10K.faa \
         -i {input} \
         -f {config[prodigal_params][output]} -p meta  > {output}
        """

rule RPSblast:
    input:"annotations/proteins/megahit_c10K.gff"
    output:"annotations/cog-annotations/megahit_c10K.out"
    shell:
        """
        set +u;source activate concoct_env;set -u
        export COGSDB_DIR={config[RPSblast_params][cogsdir]}
        {config[paths][CONCOCT]}/scripts/RPSBLAST.sh -f {config[RPSblast_params][f]} -p -c {config[RPSblast_params][c]} -r {config[RPSblast_params][r]}
        mkdir -p {config[paths][concoct_run]}/{config[RPSblast_params][outdir]}
        mv megahit_c10K.out {config[RPSblast_params][outdir]}/
        """

rule COGfilter:
    input:"annotations/cog-annotations/megahit_c10K.out","concoct-output/clustering_gt1000.csv"
    output:"evaluation-output/clustering_gt1000_scg.tab"
    shell:
        """
        set +u;source activate bcbio;set -u
        cd {config[paths][concoct_run]}
        {config[paths][CONCOCT]}/scripts/COG_table.py -b {config[COGfilter_params][b]} \
        -m {config[paths][CONCOCT]}/{config[COGfilter_params][m]} \
        -c {config[COGfilter_params][c]} \
        --cdd_cog_file {config[paths][CONCOCT]}/{config[COGfilter_params][cdd2cog]} > {output}
        """

#rule COGeval:
#    input:"evaluation-output/clustering_gt1000_scg.tab"
#    output:"evaluation-output/clustering_gt1000_scg.pdf"
#    shell:
#        """
#        set +u;source activate concoct_env;set -u
#        Rscript {config[paths][CONCOCT]}/scripts/COGPlot.R -s {input} -o {output}
#        """

#rule linkage:
#    input:"concoct-input/concoct_linkage.tsv"
#    output:"concoct-output/clustering_gt1000_l.csv"
#    shell:
#        """
#        {config[paths][CONCOCT]}/scripts/ClusterLinkNOverlap.pl --cfile={config[paths][concoct_run]}/{config[linkage_params][c]} --lfile={config[paths][concoct_run]}/{input} --covfile={config[linkage_params][cov]} --ofile={output}
#        """

rule speciesProt:
    input:"evaluation-output/clustering_gt1000_scg.tab"
    output: dynamic("carvemeOut/{species}.txt")
    shell:
        """
        cd {config[paths][concoct_run]}
        #touch {output}
        mkdir -p {config[speciesProt_params][dir]}
        cp {input} {config[paths][concoct_run]}/{config[speciesProt_params][dir]}
        cd {config[speciesProt_params][dir]}
        sed -i '1d' {config[speciesProt_params][infile]} #removes first row
        awk '{{print $2}}' {config[speciesProt_params][infile]} > allspecies.txt #extracts node information
        sed '/^>/ s/ .*//' {config[speciesProt_params][metaFASTA]} > {config[speciesProt_params][metaFASTAcleanID]} #removes annotation to protein ID
        Rscript {config[speciesProt_params][scriptdir]}multiFASTA2speciesFASTA.R
        sed -i 's/"//g' species*
        sed -i '/k99/s/^/>/' species*
        sed -i 's/{config[speciesProt_params][tab]}/{config[speciesProt_params][newline]}/' species*
        cd {config[paths][concoct_run]}
        mkdir -p {config[carveme_params][dir]}
        cd {config[carveme_params][dir]}
        cp {config[paths][concoct_run]}/{config[speciesProt_params][dir]}/species* {config[paths][concoct_run]}/{config[carveme_params][dir]}
        find . -name "species*" -size -{config[carveme_params][cutoff]} -delete #delete files with little information, these cause trouble
        cd {config[paths][concoct_run]}
        """

rule carveme:
    input: "carvemeOut/{species}.txt"
    output: "carvemeOut/{species}.xml"
    shell:
        """
        set +u;source activate concoct_env;set -u
        #echo {input}
        #echo {output}
        carve {config[carveme_params][dir]}/$(basename {input})
        #mv {output}
        #touch {output}
        """
rule memote:
    input: "carvemeOut/{species}.xml"
    output: "carvemeOut/{species}.xml.html"
    shell:
        """
        set +u;source activate concoct_env;set -u
        #cd {config[paths][concoct_run]}/{config[carveme_params][dir]}
        memote report snapshot --filename "{input}.html" {input}
        #cd {config[paths][concoct_run]}
        """

rule fetchMG:
    input: dna="megahit_c10K.fa" ,protein="annotations/proteins/megahit_c10K.faa"
    output: "fetchMG/output"
    shell:
        """
        cp -r {config[fetchMG_params][dir]} {config[paths][concoct_run]}
        cd {config[paths][concoct_run]}/{config[fetchMG_params][foldername]}
        ./fetchMG.pl -m extraction -x bin ../{input.protein} -d ../{input.dna}
        cd {config[paths][concoct_run]}
        """
rule specI:
    input: dna="megahit_c10K.fa" ,protein="annotations/proteins/megahit_c10K.faa"
    output: "specI/output"
    shell:
        """
        echo {input}
        """

#TRASH

#rule velvet:
#    input:R1="All_R1.fa", R2="All_R2.fa"
#    output:"velveth_k71/contigs.fa"
#    shell:
#        """
#        #NOTE: velvet is not the best tool to use, replace with megahit or similar tool
#        velveth {config[velvet_params][output_dir]} {config[velvet_params][hash_length]} -{config[velvet_params][file_format]} -{config[velvet_params][read_type]} \
#        -{config[velvet_params][read_type]} {input.R1} {input.R2}; velvetg {config[velvet_params][output_dir]} -ins_length {config[velvet_params][ins_length]} \
#        -exp_cov {config[velvet_params][exp_cov]} -cov_cutoff {config[velvet_params][cov_cutoff]};mkdir -p {config[velvet_params][clean_output]};
#        """

#rule cleanupvelvet:
#    input:"velveth_k71/contigs.fa"
#    output:"contigs/velvet_71.fa"
#    shell:"cp {input} {output}; rm All_R1.fa; rm All_R2.fa"

#rule bowtieVelvet:
#    input:"contigs/velvet_71_c10K.fa"
#    output:"map/Sample118_s1e5_R1.fa"
#    shell:
#        """
#        export MRKDUP={config[bowtie_params][MRKDUP_jardir]};
#        set +u;source activate concoct_env;set -u;
#        bowtie2-build {input} {input};
#        for f in map/Sample118_s1e5_R1.fa; do
#            mkdir -p map/$(basename $f);
#            cd map/$(basename $f);
#            bash {config[bowtie_params][MRKDUP_shelldir]} -ct 1 -p '-f' $f $(echo $f | sed s/R1/R2/) pair {config[paths][concoct_run]}/{input} asm bowtie2;
#            cd ../..;
#        done
#        """
