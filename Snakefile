configfile: "config.yaml"

rule all:
    input:
        "evaluation-output/ClusterPlot.pdf","evaluation-output/clustering_gt1000_scg.pdf","concoct-output/clustering_gt1000_l.csv"
    shell:"snakemake --dag | dot -Tpng > pipemap.png"

rule mergeR1:
    output:"All_R1.fa"
    shell:"cat {config[paths][CONCOCT_TEST]}/reads/Sample*_R1.fa > {output}"

rule mergeR2:
    output:"All_R2.fa"
    shell:"cat {config[paths][CONCOCT_TEST]}/reads/Sample*_R2.fa > {output}"

rule velvet:
    input:R1="All_R1.fa", R2="All_R2.fa"
    output:"velveth_k71/contigs.fa"
    shell:"velveth velveth_k71 71 -fasta -shortPaired -separate {input.R1} {input.R2}; velvetg velveth_k71 -ins_length 400 -exp_cov auto -cov_cutoff auto; mkdir contigs"

rule cleanupvelvet:
    input:"velveth_k71/contigs.fa"
    output:"contigs/velvet_71.fa"
    shell:"cp {input} {output}; rm All_R1.fa; rm All_R2.fa"

rule cutcontigs:
    input:"contigs/velvet_71.fa"
    output:"contigs/velvet_71_c10K.fa"
    shell:"set +u;source activate concoct_env;set -u;python {config[paths][CONCOCT]}/scripts/cut_up_fasta.py -c 10000 -o 0 -m {input} > {output}"

rule bowtie:
    input:"contigs/velvet_71_c10K.fa"
    output:"map/Sample343_s1e5_R1.fa/bowtie2/asm_pair-smds.bam"
    shell:
        """
        export MRKDUP=/c3se/NOBACKUP/groups/c3-c3se605-17-8/projects_francisco/.conda/envs/concoct_env/share/picard-2.18.4-0/MarkDuplicates.jar;
        set +u;source activate concoct_env;set -u;
        bowtie2-build {input} {input};
        for f in {config[paths][CONCOCT_TEST]}/reads/*_R1.fa; do
            mkdir -p map/$(basename $f);
            cd map/$(basename $f);
            bash {config[paths][CONCOCT]}/scripts/map-bowtie2-markduplicates.sh -ct 1 -p '-f' $f $(echo $f | sed s/R1/R2/) pair {config[paths][CONCOCT_EXAMPLE]}/{input} asm bowtie2;
            cd ../..;
        done
        """

rule covtable:
    input:"map/Sample343_s1e5_R1.fa/bowtie2/asm_pair-smds.bam"
    output:"concoct-input/concoct_inputtable.tsv"
    shell:
        """
        set +u;source activate concoct_env;set -u;
        cd {config[paths][CONCOCT_EXAMPLE]}/map
        python {config[paths][CONCOCT]}/scripts/gen_input_table.py --isbedfiles \
            --samplenames <(for s in Sample*; do echo $s | cut -d'_' -f1; done) \
            ../contigs/velvet_71_c10K.fa */bowtie2/asm_pair-smds.coverage > concoct_inputtable.tsv;
        mkdir -p {config[paths][CONCOCT_EXAMPLE]}/concoct-input;
        mv concoct_inputtable.tsv {config[paths][CONCOCT_EXAMPLE]}/concoct-input/;
        """

rule linktable:
    input:"contigs/velvet_71_c10K.fa"
    output:"concoct-input/concoct_linkage.tsv"
    shell:
        """
        set +u;source activate concoct_env;set -u;
        cd {config[paths][CONCOCT_EXAMPLE]}/map;
        python {config[paths][CONCOCT]}/scripts/bam_to_linkage.py -m 8 \
        --regionlength 500 --fullsearch \
        --samplenames <(for s in Sample*; do echo $s | cut -d'_' -f1; done) \
        ../{input} Sample*/bowtie2/asm_pair-smds.bam > concoct_linkage.tsv;
        mv concoct_linkage.tsv {config[paths][CONCOCT_EXAMPLE]}/concoct-input/;
        """

rule concoct:
    input:"concoct-input/concoct_inputtable.tsv"
    output:"concoct-output/clustering_gt1000.csv"
    shell:
        """
        set +u;source activate concoct_env;set -u;
        cut -f1,3- {input} > concoct-input/concoct_inputtableR.tsv;
        concoct -c 40 --coverage_file concoct-input/concoct_inputtableR.tsv --composition_file contigs/velvet_71_c10K.fa -b concoct-output/
        """

rule evaluate:
    input:"concoct-output/clustering_gt1000.csv"
    output:"evaluation-output/ClusterPlot.pdf"
    shell:
        """
        set +u;source activate concoct_env;set -u
        mkdir -p evaluation-output;
        Rscript {config[paths][CONCOCT]}/scripts/ClusterPlot.R -c {input} -p concoct-output/PCA_transformed_data_gt1000.csv -m concoct-output/pca_means_gt1000.csv -r concoct-output/pca_variances_gt1000_dim -l -o {output};
        cp {config[paths][CONCOCT_TEST]}/evaluation-output/clustering_gt1000_s.csv evaluation-output/;
        {config[paths][CONCOCT]}/scripts/Validate.pl --cfile=concoct-output/clustering_gt1000.csv --sfile=evaluation-output/clustering_gt1000_s.csv --ofile=evaluation-output/clustering_gt1000_conf.csv --ffile=contigs/velvet_71_c10K.fa;
#heatmapdont work        Rscript {config[paths][CONCOCT]}/scripts/ConfPlot.R  -c evaluation-output/clustering_gt1000_conf.csv -o  evaluation-output/clustering_gt1000_conf.pdf;
        """

rule prodigal:
    input:"contigs/velvet_71_c10K.fa"
    output:"annotations/proteins/velvet_71_c10K.gff"
    shell:
        """
        set +u;source activate concoct_env;set -u
        cd {config[paths][CONCOCT_EXAMPLE]}
        mkdir -p {config[paths][CONCOCT_EXAMPLE]}/annotations/proteins
        prodigal -a annotations/proteins/velvet_71_c10K.faa \
         -i {input} \
         -f gff -p meta  > {output}
        """

rule RPSblast:
    input:"annotations/proteins/velvet_71_c10K.gff"
    output:"annotations/cog-annotations/velvet_71_c10K.out"
    shell:
        """
        set +u;source activate concoct_env;set -u
        export COGSDB_DIR=/c3se/NOBACKUP/groups/c3-c3se605-17-8/projects_francisco/binning/snakemake-concot/cog_le
        {config[paths][CONCOCT]}/scripts/RPSBLAST.sh -f annotations/proteins/velvet_71_c10K.faa -p -c 8 -r 1
        mkdir -p {config[paths][CONCOCT_EXAMPLE]}/annotations/cog-annotations
        mv velvet_71_c10K.out annotations/cog-annotations/
        """

rule COGfilter:
    input:"annotations/cog-annotations/velvet_71_c10K.out","concoct-output/clustering_gt1000.csv"
    output:"evaluation-output/clustering_gt1000_scg.tab"
    shell:
        """
        set +u;source activate bcbio;set -u
        cd {config[paths][CONCOCT_EXAMPLE]}
        {config[paths][CONCOCT]}/scripts/COG_table.py -b annotations/cog-annotations/velvet_71_c10K.out \
        -m {config[paths][CONCOCT]}/scgs/scg_cogs_min0.97_max1.03_unique_genera.txt \
        -c concoct-output/clustering_gt1000.csv \
        --cdd_cog_file {config[paths][CONCOCT]}/scgs/cdd_to_cog.tsv > {output}
        """

rule COGeval:
    input:"evaluation-output/clustering_gt1000_scg.tab"
    output:"evaluation-output/clustering_gt1000_scg.pdf"
    shell:
        """
        set +u;source activate concoct_env;set -u
        Rscript {config[paths][CONCOCT]}/scripts/COGPlot.R -s {input} -o {output}
        """

rule linkage:
    input:"concoct-input/concoct_linkage.tsv"
    output:"concoct-output/clustering_gt1000_l.csv"
    shell:
        """
        {config[paths][CONCOCT]}/scripts/ClusterLinkNOverlap.pl --cfile={config[paths][CONCOCT_EXAMPLE]}/concoct-output/clustering_gt1000.csv --lfile={config[paths][CONCOCT_EXAMPLE]}/{input} --covfile=concoct-input/concoct_inputtableR.tsv --ofile={output}
        """

