configfile: "config.yaml"

import os
import glob

R1 = sorted([os.path.splitext(val)[0] for val in (glob.glob('/c3se/NOBACKUP/groups/c3-c3se605-17-8/projects_francisco/binning/perSamplePipe/reads/*_1.fastq'))]) #These commands dont recognize the config.yaml file
R2 = sorted([os.path.splitext(val)[0] for val in glob.glob('/c3se/NOBACKUP/groups/c3-c3se605-17-8/projects_francisco/binning/perSamplePipe/reads/*_2.fastq')]) #Change the filepath to wherever raw reads are stored
namesR1 = [os.path.basename(val) for val in R1]
names = [x[:-2] for x in namesR1]

rule all:
    input:
        expand("concoct_input/{names}/concoct_inputtableR.tsv", names = names)
        #expand("{R1files}.fastq", R1files= R1),expand("{R2files}.fastq", R2files= R2),
        #dynamic("carvemeOut/{speciesProt}.xml.html"),
        #"requiredDummyFile"
    output:
        "pipemap.png"
    shell:
        """
        snakemake --dag | dot -Tpng > {output}
        """

rule megahit:
    input:
        R1=config["paths"]["raw_reads"]+"/{names}_1.fastq",
        R2=config["paths"]["raw_reads"]+"/{names}_2.fastq"
    output:
        "megahitAssembly/{names}/final.contigs.fa"
    shell:
        """
        rm -r $(dirname {output})
        megahit -t {config[megahit_params][threads]} --tmp-dir $TMPDIR --presets meta-large --verbose -1 {input.R1} -2 {input.R2} --continue -o {config[paths][concoct_run]}/$(dirname {output})
        """

rule cutcontigs:
    input:
        "megahitAssembly/{names}/final.contigs.fa"
    output:
        "contigs/{names}/megahit_c10K.fa"
    shell:
        """
        set +u;source activate concoct_env;set -u;
        python {config[cutcontigs_params][script_dir]} -c 10000 -o 0 -m {input} > $(basename {output});
        mv $(basename {output}) {config[paths][concoct_run]}/$(dirname {output});
        """

rule kallistoBuild:
    input:
        "contigs/{names}/megahit_c10K.fa"
    output:
        "quantification/kallistoIndices/{names}.kaix"
    shell:
        """
        set +u;source activate concoct_env;set -u;
        kallisto index {input} -i {output}
        """

rule kallistoQuant:
    input:
        R1=config["paths"]["raw_reads"]+"/{names}_1.fastq",
        R2=config["paths"]["raw_reads"]+"/{names}_2.fastq",
        index="quantification/kallistoIndices/{names}.kaix"
    output:
        "quantification/{names}/abundance.tsv.gz"
    shell:
        """
        set +u;source activate concoct_env;set -u;
        kallisto quant --threads {config[kallisto_params][threads]} --plaintext -i {input.index} -o $(dirname {output}) {input.R1} {input.R2}
        gzip $(dirname {output})/abundance.tsv
        """

rule kallisto2concoctTable:
    input:
        "quantification/{names}/abundance.tsv.gz"
    output:
        "concoct_input/{names}/concoct_inputtableR.tsv"
    shell:
        """
        set +u;source activate concoct_env;set -u;
        python {config[kallisto_params][script]} \
            --samplenames <(for s in $(dirname {input}); do echo $s; done) \
                {input} > {output}
        """

rule concoct:
    input:
        table="concoct_input/{names}/concoct_inputtableR.tsv",
        comp="contigs/{names}/megahit_c10K.fa"
    output:
        "concoct_output/{names}/clustering_gt1000.csv"
    shell:
        """
        set +u;source activate concoct_env;set -u;
        cd {config[paths][concoct_run]};
        concoct --coverage_file {input.table} --composition_file {input.comp} -b $(dirname {output}) -t {config[concoct_params][threads]} -c {config[concoct_params][clusters]}
        cut -d"," -f2 $(dirname {output})/clustering_gt1000.csv | sort | uniq -c | wc > $(dirname {output})/binfo.txt
        """
        #cut {config[concoct_params][cutoptions]} {input} > {config[concoct_params][covfile]};

rule mergeClustering:
    input:
        "concoct_output/{names}/clustering_gt1000.csv"
    output:
        "concoct_output/{names}/clustering_merged.csv"
    shell:
        """
        set +u;source activate concoct_env;set -u;
        merge_cutup_clustering.py {input} > {output}
        """

rule extractBins:
    input:
        clustering="concoct_output/{names}/clustering_merged.csv",
        OGcontigs="megahitAssembly/{names}/final.contigs.fa"
    output:
        "bins/{names}"
    shell:
        """
        set +u;source activate concoct_env;set -u;
        mkdir -p {output}
        extract_fasta_bins.py {input.OGcontigs} {input.clustering} --output_path {output}
        """

rule checkM:
    input:
        "bins/{names}"
    output:
        "checkm_out/{names}"
    shell:
        """
        set +u;source activate checkm_env;set -u;
        checkm lineage_wf -x fa -t {config[checkM_params][threads]} --pplacer_threads {config[checkM_params][threads]} {input} {output}
        """

#bin filter step

#kraken2 for classification

#bracken for abundance

rule carveme:
    input: "carvemeOut/{speciesProt}.txt"
    output: "carvemeOut/{speciesProt}.xml"
    shell:
        """
        set +u;source activate concoct_env;set -u
        carve {config[carveme_params][dir]}/$(basename {input})
        """
rule memote:
    input: "carvemeOut/{speciesProt}.xml"
    output: "carvemeOut/{speciesProt}.xml.html"
    shell:
        """
        set +u;source activate concoct_env;set -u
        memote report snapshot --filename "{input}.html" {input} #generate .html report
        memote run {input} #generate quick printout of model summary
        """

                            #BELOW ARE RULES USEFUL FOR EXTRACTING BINS USING MY OWN R SCRIPT, NO LONGER USED

#rule parseFASTA:
#    input:"evaluation-output/clustering_gt1000_scg.tab"
#    output: "speciesProt"#dynamic("carvemeOut/{speciesProt}.txt"),dynamic("speciesDNA/{speciesDNA}.txt")
#    shell:
#        """
#        cd {config[paths][concoct_run]}
#        mkdir -p {config[speciesProt_params][dir]}
#        cp {input} {config[paths][concoct_run]}/{config[speciesProt_params][dir]}
#        cd {config[speciesProt_params][dir]}
#        sed -i '1d' {config[speciesProt_params][infile]} #removes first row
#        awk '{{print $2}}' {config[speciesProt_params][infile]} > allspecies.txt #extracts node information
#        sed '/^>/ s/ .*//' {config[speciesDNA_params][FASTA]} > {config[speciesDNA_params][FASTAcleanID]} #removes annotation to gene ID
#        sed '/^>/ s/ .*//' {config[speciesProt_params][pFASTA]} > {config[speciesProt_params][pFASTAcleanID]} #removes annotation to protein ID
#        Rscript {config[speciesProt_params][scriptdir]}multiFASTA2speciesFASTA.R
#        sed -i 's/"//g' species*
#        sed -i '/k99/s/^/>/' species*
#        sed -i 's/{config[speciesProt_params][tab]}/{config[speciesProt_params][newline]}/' species*
#        sed -i '1d' speciesDNA* #remove first row
#        sed -i '1d' speciesDNA* #remove second row
#        cd {config[paths][concoct_run]}
#        """

#rule pseudoProt:
#    input: "speciesProt"
#    output: dynamic("carvemeOut/{speciesProt}.txt")
#    shell:
#        """
#        cd {config[paths][concoct_run]}
#        mkdir -p {config[carveme_params][dir]}
#        cd {config[carveme_params][dir]}
#        cp {config[paths][concoct_run]}/{config[speciesProt_params][dir]}/speciesProt*.txt {config[paths][concoct_run]}/{config[carveme_params][dir]}
#        find . -name "species*.txt" -size -{config[carveme_params][cutoff]} -delete #delete files with little information, these cause trouble
#        cd {config[paths][concoct_run]}
#        """

#rule pseudoDNA:
#    input: "speciesProt"
#    output: dynamic("speciesDNA/{speciesDNA}.fa")
#    shell:
#        """
#        cd {config[paths][concoct_run]}
#        mkdir -p {config[speciesDNA_params][dir]}
#        cp {config[speciesProt_params][dir]}/speciesDNA*.txt {config[paths][concoct_run]}/{config[speciesDNA_params][dir]}
#        cp {config[speciesProt_params][dir]}/cleanID.fa {config[paths][concoct_run]}/{config[speciesDNA_params][dir]}
#        cd {config[speciesDNA_params][dir]}
#        find . -name "species*.txt" -size -{config[speciesDNA_params][cutoff]} -delete #delete files with little information, these cause trouble
#        for file in *.txt; do
#            mv "$file" "$(basename "$file" .txt).fa"
#        done
#        cd {config[paths][concoct_run]}
#        """

#rule requiredRule:
#    input: dynamic("speciesDNA/{speciesDNA}.txt")
#    output: "requiredDummyFile"
#    shell: "cd {config[paths][concoct_run]} ;touch {output}"



                ##BELOW ARE THE RULES FOR GENERATING BOWTIE INDEX, BOWTIE ALLIGNMENT, AND COVTABLE

#rule bowtieBuild:
#    input:"contigs/{names}/megahit_c10K.fa"
#    output:"contigs/{names}/buildDummy.txt"
#    shell:
#        """
#        set +u;source activate concoct_env;set -u;
#        cd $(dirname {input});
#        bowtie2-build $(basename {input}) $(basename {input});
#        cd {config[paths][concoct_run]};
#        touch {output}
#        """

#rule bowtie:
#    input:
#        assembly="contigs/{names}/megahit_c10K.fa",
#        reads=config["paths"]["raw_reads"]+"/{names}_1.fastq",
#        index= "contigs/{names}/buildDummy.txt"
#    output:"map/{names}"
#    shell:
#        """
#        set +u;source activate concoct_env;set -u;
#        export MRKDUP={config[bowtie_params][MRKDUP_jardir]};
#        DIRECTORY={output}
#        if [ ! -d "$DIRECTORY" ]; then
#            mkdir -p {output}
#            cd {output};
#            bash {config[bowtie_params][MRKDUP_shelldir]} -c -t {config[bowtie_params][threads]} -p '-q' {input.reads} $(echo {input.reads} | sed s/_1.fastq/_2.fastq/) pair {config[paths][concoct_run]}/{input.assembly} asm bowtie2;
#            cd {config[paths][concoct_run]};
#        fi
#        """

#rule covtable:
#    input: "map/{names}"
#    output:"concoct-input/{names}/concoct_inputtable.tsv"
#    shell:
#        """
#        set +u;source activate concoct_env;set -u;
#        cd {config[paths][concoct_run]}/{input}
#        python {config[paths][CONCOCT]}/scripts/gen_input_table.py --isbedfiles \
#            --samplenames <(for s in ERR*; do echo $s | cut -d'_' -f1; done) \
#            ../contigs/{names}/megahit_c10K.fa */bowtie2/asm_pair-smds.coverage > concoct_inputtable.tsv;
#        mkdir -p {config[paths][concoct_run]}/concoct_input/{names};
#        mv concoct_inputtable.tsv {config[paths][concoct_run]}/concoct_input/{names};
#        """
