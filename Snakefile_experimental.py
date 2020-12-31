rule metaspades: 
    input:
        R1 = rules.qfilter.output.R1, 
        R2 = rules.qfilter.output.R2
    output:
        config["path"]["root"]+"/"+config["folder"]["assemblies"]+"/{IDs}/contigs.fasta.gz"
    benchmark:
        config["path"]["root"]+"/"+"benchmarks/{IDs}.metaspades.benchmark.txt"
    shell:
        """
        set +u;source activate {config[envs][metabagpipes]};set -u;
        cp {input.R1} {input.R2} $TMPDIR    
        cd $TMPDIR
        metaspades.py --only-assembler -1 $(basename {input.R1}) -2 $(basename {input.R2}) -t {config[cores][metaspades]} -o .
        gzip contigs.fasta
        mkdir -p $(dirname {output})
        rm $(basename {input.R1}) $(basename {input.R2})
        mv -v contigs.fasta.gz spades.log $(dirname {output})
        """

rule metabatMultiSample:
    input:
        assembly=config["path"]["root"]+"/"+config["folder"]["assemblies"]+"/{IDs}/contigs.fasta.gz",
        reads=config["path"]["root"]+"/"+config["folder"]["qfiltered"]
    output:
        directory(config["path"]["root"]+"/"+config["folder"]["metabat"]+"/{IDs}/{IDs}.metabat-bins")
    benchmark:
        config["path"]["root"]+"/"+"benchmarks/{IDs}.metabat.benchmark.txt"
    shell:
        """
        set +u;source activate {config[envs][metabagpipes]};set -u;
        mkdir -p $(dirname $(dirname {output}))
        mkdir -p $(dirname {output})
        cp {input.assembly} $TMPDIR
        cd $TMPDIR
        mv $(basename {input.assembly}) $(basename $(dirname {input.assembly})).gz
        gunzip $(basename $(dirname {input.assembly})).gz
        bwa index $(basename $(dirname {input.assembly}))

        for sample in {input.reads}/*;do
            echo "Mapping sample $sample ... "
            ID=$(basename $sample);
            bwa mem -t {config[cores][metabat]} $(basename $(dirname {input.assembly})) $sample/*_1.fastq.gz $sample/*_2.fastq.gz > $ID.sam
            samtools view -@ {config[cores][metabat]} -Sb $ID.sam > $ID.bam
            samtools sort -@ {config[cores][metabat]} $ID.bam $ID.sort
            rm $ID.bam $ID.sam
            echo "Done mapping sample $sample !"
            echo "Creating depth file for sample $sample ... "
            jgi_summarize_bam_contig_depths --outputDepth depth.txt $ID.sort.bam
            echo "Done creating depth file for sample $sample !"
            rm $ID.sort.bam
            paste $sample.depth.txt
        done

        runMetaBat.sh $(basename $(dirname {input.assembly})) *.sort.bam
        mv *.txt *.tab $(basename {output}) $(dirname {output})
        """

rule maxbinMultiSample:
    input:
        assembly=config["path"]["root"]+"/"+config["folder"]["assemblies"]+"/{IDs}/contigs.fasta.gz",
        reads=config["path"]["root"]+"/"+config["folder"]["qfiltered"]
    output:
        directory(config["path"]["root"]+"/"+config["folder"]["maxbin"]+"/{IDs}/{IDs}.maxbin-bins")
    benchmark:
        config["path"]["root"]+"/"+"benchmarks/{IDs}.maxbin.benchmark.txt"
    shell:
        """
        set +u;source activate {config[envs][metabagpipes]};set -u;

        mkdir -p $(dirname $(dirname {output}))
        mkdir -p $(dirname {output})
        cp {input.assembly} $TMPDIR
        cd $TMPDIR

        focal=$(basename $(dirname {input.assembly}))
        gunzip contigs.fasta.gz

        echo "Creating kallisto index for focal sample $focal ... "
        #kallisto index contigs.fasta -i $focal.kaix
        cp /home/zorrilla/workspace/straya/test/*.kaix .
        echo "Done creating kallisto index!"

        echo "Begin cross mapping samples ... "
        for folder in {input.reads}/*;do
            var=$(basename $folder)
            echo "Mapping sample $var to focal sample $focal using kallisto quant ... "
            kallisto quant --threads {config[cores][kallisto]} --plaintext -i $focal.kaix -o . $folder/*_1.fastq.gz $folder/*_2.fastq.gz;
            #tail -n +2 abundance.tsv > $(basename $folder)_abundance.tsv
            #rm abundance.tsv
            echo "Done mapping sample $var to focal sample!"
        done
        echo "Done cross mapping all samples! "

        find . -name "*_abundance.tsv" > abund_list.txt

        echo "Begin running maxbin2 algorithm ... "
        run_MaxBin.pl -contig contigs.fasta -out $focal -abund_list abund_list.txt -thread {config[cores][maxbin]}
        echo "Done running maxbin2!"

        rm contigs.fasta
        mkdir $(basename {output})
        mv *.fasta $(basename {output})
        mv $(basename {output}) $(dirname {output})
        """ 

rule mOTUs2classifyGenomes:
    input:
        bins = f'{config["path"]["root"]}/{config["folder"]["reassembled"]}/{{IDs}}/reassembled_bins',
        script = f'{config["path"]["root"]}/{config["folder"]["scripts"]}/classify-genomes'
    output:
        #directory(f'{config["path"]["root"]}/{config["folder"]["classification"]}/{{IDs}}')
    benchmark:
        f'{config["path"]["root"]}/benchmarks/{{IDs}}.classify-genomes.benchmark.txt'
    shell:
        """
        set +u;source activate {config[envs][metabagpipes]};set -u;
        mkdir -p {output}
        cd $SCRATCHDIR
        cp -r {input.script}/* {input.bins}/* .

        echo "Begin classifying bins ... "
        for bin in *.fa; do
            echo -e "\nClassifying $bin ... "
            $PWD/classify-genomes $bin -t {config[cores][classify]} -o $(echo $bin|sed 's/.fa/.taxonomy/')
            cp *.taxonomy {output}
            rm *.taxonomy
            rm $bin 
        done
        echo "Done classifying bins. "
        """

rule taxonomyVis:
    input: 
        f'{config["path"]["root"]}/{config["folder"]["classification"]}'
    output: 
        text = f'{config["path"]["root"]}/{config["folder"]["stats"]}/classification.stats',
        plot = f'{config["path"]["root"]}/{config["folder"]["stats"]}/taxonomyVis.pdf'
    message:
        """
        mOTUs2 taxonomy visualization.
        Generate bar plot with most common taxa (n>15) and density plots with mapping statistics.
        """
    shell:
        """
        set +u;source activate {config[envs][metabagpipes]};set -u;
        cd {input}

        echo -e "\nBegin reading classification result files ... \n"
        for folder in */;do 

            for file in $folder*.taxonomy;do

                # Define sample ID to append to start of each bin name in summary file
                sample=$(echo $folder|sed 's|/||')

                # Define bin name with sample ID, shorten metaWRAP naming scheme (orig/permissive/strict)
                fasta=$(echo $file | sed 's|^.*/||' | sed 's/.taxonomy//g' | sed 's/orig/o/g' | sed 's/permissive/p/g' | sed 's/strict/s/g' | sed "s/^/$sample./g");

                # Extract NCBI ID 
                NCBI=$(less $file | grep NCBI | cut -d ' ' -f4);

                # Extract consensus taxonomy
                tax=$(less $file | grep tax | sed 's/Consensus taxonomy: //g');

                # Extract consensus motus
                motu=$(less $file | grep mOTUs | sed 's/Consensus mOTUs: //g');

                # Extract number of detected genes
                detect=$(less $file | grep detected | sed 's/Number of detected genes: //g');

                # Extract percentage of agreeing genes
                percent=$(less $file | grep agreeing | sed 's/Percentage of agreeing genes: //g' | sed 's/%//g');

                # Extract number of mapped genes
                map=$(less $file | grep mapped | sed 's/Number of mapped genes: //g');
                
                # Extract COG IDs, need to use set +e;...;set -e to avoid erroring out when reading .taxonomy result file for bin with no taxonomic annotation
                set +e
                cog=$(less $file | grep COG | cut -d$'\t' -f1 | tr '\n' ',' | sed 's/,$//g');
                set -e
                
                # Display and store extracted results
                echo -e "$fasta \t $NCBI \t $tax \t $motu \t $detect \t $map \t $percent \t $cog"
                echo -e "$fasta \t $NCBI \t $tax \t $motu \t $detect \t $map \t $percent \t $cog" >> classification.stats;
            
            done;
        
        done

        echo -e "\nDone generating classification.stats summary file, moving to stats/ directory and running taxonomyVis.R script ... "
        mv classification.stats {config[path][root]}/{config[folder][stats]}
        cd {config[path][root]}/{config[folder][stats]}

        Rscript {config[path][root]}/{config[folder][scripts]}/{config[scripts][taxonomyVis]}
        rm Rplots.pdf # Delete redundant pdf file
        echo "Done. "
        """


rule parseTaxAb:
    input:
        taxonomy = rules.taxonomyVis.output.text ,
        abundance = f'{config["path"]["root"]}/{config["folder"]["abundance"]}'
    output:
        directory(f'{config["path"]["root"]}/MAG.table')
    message:
        """
        Parses an abundance table with MAG taxonomy for rows and samples for columns.
        Note: parseTaxAb should only be run after the classifyGenomes, taxonomyVis, and abundance rules.
        """
    shell:
        """
        set +u;source activate {config[envs][metabagpipes]};set -u
        cd {input.abundance}

        for folder in */;do

            # Define sample ID
            sample=$(echo $folder|sed 's|/||g')
            
            # Same as in taxonomyVis rule, modify bin names by adding sample ID and shortening metaWRAP naming scheme (orig/permissive/strict)
            paste $sample/$sample.abund | sed 's/orig/o/g' | sed 's/permissive/p/g' | sed 's/strict/s/g' | sed "s/^/$sample./g" >> abundance.stats
       
        done

        mv abundance.stats {config[path][root]}/{config[folder][stats]}
        cd {config[path][root]}/{config[folder][stats]}

        """
