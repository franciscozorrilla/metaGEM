configfile: "config.yaml"

import os
import glob

IDs = sorted([os.path.splitext(val)[0] for val in (glob.glob('raw/*'))])
IDs = [os.path.basename(val) for val in IDs] #grab just sample ID

#Make sure that final_bins/ folder contains all bins in single folder for binIDs wildcard to work. Use moveBins rule or perform manually.
binIDs = sorted([os.path.splitext(val)[0] for val in (glob.glob('final_bins/*'))])
binIDs = [os.path.basename(val) for val in binIDs] #grab just bin ID

rule all:
    input:
        expand(config["path"]["root"]+"/"+config["folder"]["SMETANA"]+"/{IDs}.tsv", IDs = IDs)
    shell:
        """
        echo {input}
        """

rule organizeData:
    message:
        "Assuming all samples are downloaded to the same directory, sorting paired end raw reads into sample specific sub folders within the 'dataset' folder"
    shell:
        """
        cd dataset
        for file in *.gz;do echo $file;done|sed 's/_.*$//g'|uniq > ID_samples.txt
        while read line;do mkdir -p $line;mv $line*.gz $line;done < ID_samples.txt
        rm ID_samples.txt
        cd config["path"]["root"]
        """

rule createFolders:
    input:
        {config["path"]["root"]}
    message:
        "Not strictly necessary, as we could create individual folders within rules, but doing this at the beginning may give the user a better idea of what output/intermediate files to expect"
    shell:
        """
        cd {input}
        paste config.yaml |cut -d':' -f2|tail -n +4|head -n 15 > folders.txt
        while read line;do mkdir -p $line;done < folders.txt
        rm folders.txt
        """      

rule metaspades: 
    input:
        R1=config["path"]["root"]+"/"+config["folder"]["data"]+"/{IDs}/{IDs}_1.fastq.gz", 
        R2=config["path"]["root"]+"/"+config["folder"]["data"]+"/{IDs}/{IDs}_2.fastq.gz" 
    output:
        config["path"]["root"]+"/"+config["folder"]["assemblies"]+"/{IDs}/contigs.fasta.gz"
    benchmark:
        config["path"]["root"]+"/"+"benchmarks/{IDs}.metaspades.benchmark.txt"
    shell:
        """
        set +u;source activate binenv;set -u;
        cp {input.R1} {input.R2} $TMPDIR    
        cd $TMPDIR
        metaspades.py -1 $(basename {input.R1}) -2 $(basename {input.R2}) -t 48 -o .
        gzip contigs.fasta
        mkdir -p $(dirname {output})
        mv -v contigs.fasta.gz spades.log $(dirname {output})
        """

rule kallisto:
    input:
        contigs=config["path"]["root"]+"/"+config["folder"]["assemblies"]+"/{IDs}/contigs.fasta.gz",
        reads=config["path"]["root"]+"/"+config["folder"]["data"]+"/"
    output:
        config["path"]["root"]+"/"+config["folder"]["concoctInput"]+"/{IDs}_concoct_inputtableR.tsv"
    benchmark:
        config["path"]["root"]+"/"+"benchmarks/{IDs}.kallisto.benchmark.txt"
    shell:
        """
        set +u;source activate binenv;set -u;
        cp {input.contigs} $TMPDIR
        cd $TMPDIR
        gunzip $(basename {input.contigs})
        cut_up_fasta.py -c 10000 -o 0 -m contigs.fasta > metaspades_c10k.fa
        kallisto index metaspades_c10k.fa -i $(basename $(dirname {input.contigs})).kaix
        for folder in {input.reads}*;do
        cp $folder/*.fastq.gz $TMPDIR;
        kallisto quant --threads 48 --plaintext -i $(basename $(dirname {input.contigs})).kaix -o . *_1.fastq.gz *_2.fastq.gz;
        gzip abundance.tsv;
        mv abundance.tsv.gz $(echo $(basename $folder)_abundance.tsv.gz);
        rm *.fastq.gz;
        done
        python {config[scripts][kallisto2concoct]} --samplenames <(for s in *abundance.tsv.gz; do echo $s|sed 's/_abundance.tsv.gz//'g; done) *abundance.tsv.gz > $(basename {output})
        mv $(basename {output}) $(dirname {output})
        """

rule concoct:
    input:
        table=config["path"]["root"]+"/"+config["folder"]["concoctInput"]+"/{IDs}_concoct_inputtableR.tsv",
        contigs=config["path"]["root"]+"/"+config["folder"]["assemblies"]+"/{IDs}/contigs.fasta.gz"
    output:
        directory(config["path"]["root"]+"/"+config["folder"]["concoctOutput"]+"/{IDs}/{IDs}.concoct-bins")
    benchmark:
        config["path"]["root"]+"/"+"benchmarks/{IDs}.concoct.benchmark.txt"
    shell:
        """
        set +u;source activate binenv;set -u;
        mkdir -p $(dirname {output})
        mkdir -p $(dirname $(dirname {output}))
        cp {input.contigs} {input.table} $TMPDIR
        cd $TMPDIR
        gunzip $(basename {input.contigs})
        cut_up_fasta.py -c 10000 -o 0 -m contigs.fasta > metaspades_c10k.fa
        concoct --coverage_file {input.table} --composition_file metaspades_c10k.fa -b $(basename $(dirname {output})) -t 48 -c 800
        merge_cutup_clustering.py $(basename $(dirname {output}))_clustering_gt1000.csv > $(basename $(dirname {output}))_clustering_merged.csv
        mkdir -p $(basename {output})
        extract_fasta_bins.py contigs.fasta $(basename $(dirname {output}))_clustering_merged.csv --output_path $(basename {output})
        checkm lineage_wf -x fa -t 48 --pplacer_threads 48 --file $(basename $(dirname {output})).tab --tab_table $(basename {output}) .
        mv $(basename {output}) *.log *.txt *.csv *.tab $(dirname {output})
        """

rule metabat:
    input:
        assembly=config["path"]["root"]+"/"+config["folder"]["assemblies"]+"/{IDs}/contigs.fasta.gz",
        R1=config["path"]["root"]+"/"+config["folder"]["data"]+"/{IDs}/{IDs}_1.fastq.gz", 
        R2=config["path"]["root"]+"/"+config["folder"]["data"]+"/{IDs}/{IDs}_2.fastq.gz" 
    output:
        directory(config["path"]["root"]+"/"+config["folder"]["metabat"]+"/{IDs}/{IDs}.metabat-bins")
    benchmark:
        config["path"]["root"]+"/"+"benchmarks/{IDs}.metabat.benchmark.txt"
    shell:
        """
        set +u;source activate binenv;set -u;
        mkdir -p $(dirname $(dirname {output}))
        mkdir -p $(dirname {output})
        cp {input.assembly} {input.R1} {input.R2} $TMPDIR
        cd $TMPDIR
        mv $(basename {input.assembly}) $(basename $(dirname {input.assembly})).gz
        gunzip $(basename $(dirname {input.assembly})).gz
        bwa index $(basename $(dirname {input.assembly}))
        bwa mem -t 48 $(basename $(dirname {input.assembly})) $(basename {input.R1}) $(basename {input.R2}) > $(basename $(dirname {input.assembly})).sam
        samtools view -@ 48 -Sb $(basename $(dirname {input.assembly})).sam > $(basename $(dirname {input.assembly})).bam
        samtools sort -@ 48 $(basename $(dirname {input.assembly})).bam > $(basename $(dirname {input.assembly})).sort
        runMetaBat.sh $(basename $(dirname {input.assembly})) $(basename $(dirname {input.assembly})).sort
        checkm lineage_wf -x fa -t 48 --pplacer_threads 48 --file $(basename $(dirname {input.assembly})).tab --tab_table $(basename $(dirname {input.assembly})).metabat-bins .
        mv *.txt *.tab $(basename {output}) $(dirname {output})
        """

rule maxbin:
    input:
        assembly=config["path"]["root"]+"/"+config["folder"]["assemblies"]+"/{IDs}/contigs.fasta.gz",
        R1=config["path"]["root"]+"/"+config["folder"]["data"]+"/{IDs}/{IDs}_1.fastq.gz", 
        R2=config["path"]["root"]+"/"+config["folder"]["data"]+"/{IDs}/{IDs}_2.fastq.gz" 
    output:
        directory(config["path"]["root"]+"/"+config["folder"]["maxbin"]+"/{IDs}/{IDs}.maxbin-bins")
    benchmark:
        config["path"]["root"]+"/"+"benchmarks/{IDs}.maxbin.benchmark.txt"
    shell:
        """
        set +u;source activate binenv;set -u;
        mkdir -p $(dirname $(dirname {output}))
        mkdir -p $(dirname {output})
        cp {input.assembly} {input.R1} {input.R2} $TMPDIR
        cd $TMPDIR
        gunzip contigs.fasta.gz
        run_MaxBin.pl -contig contigs.fasta -out $(basename $(dirname {output})) -reads *_1.fastq.gz -reads2 *_2.fastq.gz -thread 48
        rm contigs.fasta *.gz
        mkdir $(basename {output})
        mv *.fasta $(basename {output})
        checkm lineage_wf -x fasta -t 48 --pplacer_threads 48 --file $(basename $(dirname {input.assembly})).tab --tab_table $(basename {output}) .
        mv *.abund1 *.abund2 *.tab $(basename {output}) $(dirname {output})
        """

rule binRefine:
    input:
        concoct=config["path"]["root"]+"/"+config["folder"]["concoctOutput"]+"/{IDs}/{IDs}.concoct-bins",
        metabat=config["path"]["root"]+"/"+config["folder"]["metabat"]+"/{IDs}/{IDs}.metabat-bins",
        maxbin=config["path"]["root"]+"/"+config["folder"]["maxbin"]+"/{IDs}/{IDs}.maxbin-bins"
    output:
        directory(config["path"]["root"]+"/"+config["folder"]["refined"]+"/{IDs}")
    benchmark:
        config["path"]["root"]+"/"+"benchmarks/{IDs}.binRefine.benchmark.txt"
    shell:
        """
        set +u;source activate metawrap;set -u;
        cp -r {input.concoct} {input.metabat} {input.maxbin} $TMPDIR
        mkdir -p $(dirname {output})
        cd $TMPDIR
        metaWRAP bin_refinement -o . -A $(basename {input.concoct}) -B $(basename {input.metabat}) -C $(basename {input.maxbin}) -t 48 -m 1600 -c 50 -x 10
        rm -r $(basename {input.concoct}) $(basename {input.metabat}) $(basename {input.maxbin})
        mkdir -p {output}
        mv * {output}
        """

rule binReassemble:
    input:
        R1=config["path"]["root"]+"/"+config["folder"]["data"]+"/{IDs}/{IDs}_1.fastq.gz", 
        R2=config["path"]["root"]+"/"+config["folder"]["data"]+"/{IDs}/{IDs}_2.fastq.gz",
        refinedBins=config["path"]["root"]+"/"+config["folder"]["refined"]+"/{IDs}/metawrap_50_10_bins"
    output:
        directory(config["path"]["root"]+"/"+config["folder"]["reassembled"]+"/{IDs}")
    benchmark:
        config["path"]["root"]+"/"+"benchmarks/{IDs}.binReassemble.benchmark.txt"
    shell:
        """
        set +u;source activate metawrap;set -u;
        mkdir -p $(dirname {output})
        cp -r {input.refinedBins} {input.R1} {input.R2} $TMPDIR
        cd $TMPDIR
        metaWRAP reassemble_bins -o $(basename {output}) -b $(basename {input.refinedBins}) -1 $(basename {input.R1}) -2 $(basename {input.R1}) -t 48 -m 500 -c 50 -x 10
        rm -r $(basename {input.refinedBins})
        rm -r $(basename {output})/work_files
        rm *.fastq.gz 
        mv * $(dirname {output})
        """

rule classifyGenomes: 
    input:
        bins=config["path"]["root"]+"/"+config["folder"]["reassembled"]+"/{IDs}/reassembled_bins",
        script=config["path"]["root"]+"/"+config["folder"]["scripts"]+"/classify-genomes"
    output:
        directory(config["path"]["root"]+"/"+config["folder"]["classification"]+"/{IDs}")
    benchmark:
        config["path"]["root"]+"/"+"benchmarks/{IDs}.classify-genomes.benchmark.txt"
    shell:
        """
        set +u;source activate binenv;set -u;
        mkdir -p {output}
        cd $TMPDIR
        cp -r {input.script}/* {input.bins}/* .
        for bin in *.fa; do
        echo "RUNNING BIN $bin"
        $PWD/classify-genomes $bin -t 2 -o $(echo $bin|sed 's/.fa/.taxonomy/')
        echo "DONE" 
        echo " "
        cp *.taxonomy {output}
        rm *.taxonomy
        rm $bin 
        done
        """

rule abundance:
    input:
        bins=config["path"]["root"]+"/"+config["folder"]["reassembled"]+"/{IDs}/reassembled_bins",
        R1=config["path"]["root"]+"/"+config["folder"]["data"]+"/{IDs}/{IDs}_1.fastq.gz", 
        R2=config["path"]["root"]+"/"+config["folder"]["data"]+"/{IDs}/{IDs}_2.fastq.gz" 
    output:
        directory(config["path"]["root"]+"/"+config["folder"]["abundance"]+"/{IDs}")
    benchmark:
        config["path"]["root"]+"/"+"benchmarks/{IDs}.abundance.benchmark.txt"
    message:
        """
        Calculate bin abundance fraction using the following:

        binAbundanceFraction = (L * X / Y / Z) * 100

        X = # of reads mapped to bin_i by mapping reads to bin_i using samtools.
        Y = length of bin_i.
        Z = total # of reads mapped to all bins in sample_k
        L = length of reads (100 bp)
        """
    shell:
        """
        mkdir -p {output}
        cd $TMPDIR
        cp {input.bins}/* .
        cat *.fa > $(basename {output}).fa
        cp {input.R1} {input.R2} .
        echo "CREATING INDEX FOR BIN CONCATENATION AND MAPPING READS "
        bwa index $(basename {output}).fa
        bwa mem -t 16 $(basename {output}).fa $(basename {input.R1}) $(basename {input.R2}) > $(basename {output}).sam
        samtools view -@ 16 -Sb $(basename {output}).sam > $(basename {output}).bam
        samtools sort -@ 16 $(basename {output}).bam > $(basename {output}).sort
        samtools flagstat $(basename {output}).sort > map.stats
        cp map.stats {output}/$(basename {output})_map.stats
        rm $(basename {output}).fa *.sam *.bam *.sort
        echo "DONE MAPPING READS TO BIN CONCATENATION, BEGIN MAPPING READS TO EACH BIN "
        for bin in *.fa;do
            mkdir -p $(echo "$bin"| sed "s/.fa//")
            mv $bin $(echo "$bin"| sed "s/.fa//")
            cd $(echo "$bin"| sed "s/.fa//")
            echo "INDEXING AND MAPPING BIN $bin "
            bwa index $bin
            bwa mem -t 16 $bin ../$(basename {input.R1}) ../$(basename {input.R2}) > $(echo "$bin"|sed "s/.fa/.sam/")
            samtools view -@ 16 -Sb $(echo "$bin"|sed "s/.fa/.sam/") > $(echo "$bin"|sed "s/.fa/.bam/")
            samtools sort -@ 16 $(echo "$bin"|sed "s/.fa/.bam/") > $(echo "$bin"|sed "s/.fa/.sort/")
            samtools flagstat $(echo "$bin"|sed "s/.fa/.sort/") > $(echo "$bin"|sed "s/.fa/.map/")
            echo -n "Bin Length = " >> $(echo "$bin"|sed "s/.fa/.map/")
            less $bin|cut -d '_' -f4| awk -F' ' '{{print $NF}}'|sed 's/len=//'|awk '{{sum+=$NF;}}END{{print sum;}}' >> $(echo "$bin"|sed "s/.fa/.map/")
            echo "FINISHED MAPPING READS TO BIN $bin "
            paste $(echo "$bin"|sed "s/.fa/.map/")
            echo "CALCULATING ABUNDANCE FOR BIN $bin "
            echo -n "$bin"|sed "s/.fa//" >> $(echo "$bin"|sed "s/.fa/.abund/")
            echo -n $'\t' >> $(echo "$bin"|sed "s/.fa/.abund/")
            X=$(less $(echo "$bin"|sed "s/.fa/.map/")|grep "mapped ("|awk -F' ' '{{print $1}}')
            Y=$(less $(echo "$bin"|sed "s/.fa/.map/")|tail -n 1|awk -F' ' '{{print $4}}')
            Z=$(less "../map.stats"|grep "mapped ("|awk -F' ' '{{print $1}}')
            awk -v x="$X" -v y="$Y" -v z="$Z" 'BEGIN{{print (100*x/y/z) * 100}}' >> $(echo "$bin"|sed "s/.fa/.abund/")
            rm $bin
            cp $(echo "$bin"|sed "s/.fa/.map/") {output}
            mv $(echo "$bin"|sed "s/.fa/.abund/") ../
            cd ..
        done
        cat *.abund > $(basename {output}).abund
        mv $(basename {output}).abund {output}
        """

rule moveBins:
    message:
        "Moves and renames all bins from reassembled_bins/{sampleID}/ subfolders to single folder final_bins/. Need all bins in one folder to run rules with {binIDs} wildcard."
    shell:
        """
        cd {config[path][root]}
        mkdir -p final_bins
        for folder in reassembled_bins/*/;do for bin in $folder/reassembled_bins/*.fa;do mv $bin final_bins/$(echo $bin|sed 's|reassembled_bins/||g'|sed 's|/|_|g'|sed 's/strict/s/g'|sed 's/orig/o/g'|sed 's/permissive/p/g');done;done
        """

rule carveme:
    input:
        bin=config["path"]["root"]+"/"+"final_bins/{binIDs}.fa",
        media=config["dbs"]["carveme"]
    output:
        config["path"]["root"]+"/"+config["folder"]["GEMs"]+"/{binIDs}.xml"
    benchmark:
        config["path"]["root"]+"/"+"benchmarks/{binIDs}.carveme.benchmark.txt"
    shell:
        """
        set +u;source activate carvenv;set -u
        mkdir -p $(dirname {output})
        cp {input.bin} {input.media} $TMPDIR
        cd $TMPDIR
        carve -g M8 -v --mediadb $(basename {input.media}) --fbc2 --dna $(basename {input.bin}) -o $(echo $(basename {input.bin})|sed 's/.fa/.xml/g')
        [ -f *.xml ] && mv *.xml $(dirname {output})
        """

#rule carvemeSample:
#    input:
#        bins=config["path"]["root"]+"/"+config["folder"]["reassembled"]+"/{IDs}/reassembled_bins",
#        media=config["dbs"]["carveme"]
#    output:
#        directory(config["path"]["root"]+"/"+config["folder"]["GEMs"]+"/{IDs}")
#    benchmark:
#        config["path"]["root"]+"/"+"benchmarks/{IDs}.carveme.benchmark.txt"
#    shell:
#        """
#        set +u;source activate carvenv;set -u
#        mkdir -p {output}
#        cp {input.bins}/*.fa {input.media} $TMPDIR
#        cd $TMPDIR
#        for bin in *.fa;do
#        carve -g M8 -v --mediadb $(basename {input.media}) --fbc2 --dna $bin -o $(echo $bin| sed 's/.fa/.xml/g')
#        [ -f *.xml ] && mv *.xml {output}
#        rm $bin
#        done
#        """

rule organizeGEMs:
    message:
        "Organizes GEMs into sample specific subfolders. Necessary to run smetana per sample using the {IDs} wildcard."
    shell:
        """
        cd {config[path][refined]}
        for folder in */;do mkdir -p ../{config[path][GEMs]}; mv ../{config[path][GEMs]}/$(echo $folder|sed 's|/||')_*.xml ../{config[path][GEMs]}/$folder;done
        """

rule smetana:
    input:
        config["path"]["root"]+"/"+config["folder"]["GEMs"]+"/{IDs}"
    output:
        config["path"]["root"]+"/"+config["folder"]["SMETANA"]+"/{IDs}.tsv"
    benchmark:
        config["path"]["root"]+"/"+"benchmarks/{IDs}.smetana.benchmark.txt"
    shell:
        """
        set +u;source activate memotenv;set -u
        mkdir -p {config[path][root]}/{config[folder][SMETANA]}
        cp {config[dbs][carveme]} {input}/*.xml $TMPDIR
        cd $TMPDIR
        smetana -o $(basename {input}) --flavor fbc2 --mediadb media_db.tsv -m M1,M2,M3,M4,M5,M7,M8,M9,M10,M11,M13,M14,M15A,M15B,M16 --detailed --solver CPLEX -v *.xml
        cp *.tsv {config[path][root]} #safety measure for backup of results in case rule fails for some reason
        mv *.tsv $(dirname {output})
        """

rule memote:
    input:
        config["path"]["root"]+"/"+config["folder"]["GEMs"]+"/{IDs}"
    output:
        directory(config["path"]["root"]+"/"+config["folder"]["memote"]+"/{IDs}")
    benchmark:
        config["path"]["root"]+"/"+"benchmarks/{IDs}.memote.benchmark.txt"
    shell:
        """
        set +u;source activate memotenv;set -u
        mkdir -p $(dirname {output})
        mkdir -p {output}
        cp {input}/*.xml $TMPDIR
        cd $TMPDIR
        for model in *.xml;do
        memote report snapshot --filename $(echo $model|sed 's/.xml/.html/') $model
        memote run $model > $(echo $model|sed 's/.xml/-summary.txt/')
        mv *.txt *.html {output}
        rm $model
        done
        """

rule grid:
    input:
        bins=config["path"]["root"]+"/"+config["folder"]["reassembled"]+"/{IDs}/reassembled_bins",
        R1=config["path"]["root"]+"/"+config["folder"]["data"]+"/{IDs}/{IDs}_1.fastq.gz", 
        R2=config["path"]["root"]+"/"+config["folder"]["data"]+"/{IDs}/{IDs}_2.fastq.gz" 
    output:
        directory(config["path"]["root"]+"/"+config["folder"]["GRiD"]+"/{IDs}")
    benchmark:
        config["path"]["root"]+"/"+"benchmarks/{IDs}.grid.benchmark.txt"
    shell:
        """
        set +u;source activate binenv;set -u
        cp -r {input.bins} {input.R1} {input.R2} $TMPDIR
        cd $TMPDIR
        mkdir MAGdb
        update_database -d MAGdb -g $(basename {input.bins}) -p MAGdb
     	rm -r $(basename {input.bins})
        """

rule gridTest:
    input:
        db="/home/zorrilla/workspace/tutorial/test/MAGdb",
        sample="/home/zorrilla/workspace/tutorial/test/ERR260173.fastq.gz"
    output:
        "/home/zorrilla/workspace/tutorial/test/out3"
    shell:
        """
        set +u;source activate binenv;set -u
        cp -r {input.db} {input.sample} $TMPDIR
        cd $TMPDIR
        mkdir out3
        grid multiplex -r . -e fastq.gz -d MAGdb -p -c 0.2 -o out3 -n 48
        mv out3 /home/zorrilla/workspace/tutorial/test/
     	"""

rule smetanaTest:
    shell:
        """
        set +u;source activate memotenv;set -u
        cp -r /home/zorrilla/workspace/tutorial/test/smetana $TMPDIR
        cd $TMPDIR
        smetana -o TEST --flavor fbc2 --mediadb smetana/media_db.tsv -m M1,M2,M3,M4,M5,M7,M8,M9,M10,M11,M13,M14,M15A,M15B,M16 --detailed --solver CPLEX -v smetana/*.xml
        mv *.tsv /home/zorrilla/workspace/tutorial/test/smetana
     	"""

