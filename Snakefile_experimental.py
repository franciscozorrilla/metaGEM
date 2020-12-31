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

rule megahitCoassembly:
    input:
        R1 = f'/scratch/zorrilla/soil/coassembly/data/{{borkSoil}}/R1', 
        R2 = f'/scratch/zorrilla/soil/coassembly/data/{{borkSoil}}/R2'
    output:
        f'{config["path"]["root"]}/coassembly/coassemblies/{{borkSoil}}/contigs.fasta.gz'
    benchmark:
        f'{config["path"]["root"]}/benchmarks/coassembly.{{borkSoil}}.benchmark.txt'
    shell:
        """
        set +u;source activate {config[envs][metabagpipes]};set -u;
        cd $SCRATCHDIR

        echo -n "Copying qfiltered reads to $SCRATCHDIR ... "
        cp -r {input.R1} {input.R2} $SCRATCHDIR
        echo "done. "

        R1=$(ls R1/|tr '\n' ','|sed 's/,$//g')
        R2=$(ls R2/|tr '\n' ','|sed 's/,$//g')

        mv R1/* .
        mv R2/* .

        echo -n "Running megahit ... "
        megahit -t {config[cores][megahit]} \
            --presets {config[params][assemblyPreset]} \
            --min-contig-len {config[params][assemblyMin]}\
            --verbose \
            -1 $R1 \
            -2 $R2 \
            -o tmp;
        echo "done. "

        echo "Renaming assembly ... "
        mv tmp/final.contigs.fa contigs.fasta
        
        echo "Fixing contig header names: replacing spaces with hyphens ... "
        sed -i 's/ /-/g' contigs.fasta

        echo "Zipping and moving assembly ... "
        gzip contigs.fasta
        mkdir -p $(dirname {output})
        mv contigs.fasta.gz $(dirname {output})
        echo "Done. "
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


rule crossMap2:  
    input:
        contigs = f'{config["path"]["root"]}/{config["folder"]["assemblies"]}/{{focal}}/contigs.fasta.gz',
        R1 = rules.qfilter.output.R1,
        R2 = rules.qfilter.output.R2
    output:
        directory(f'{config["path"]["root"]}/{config["folder"]["crossMap"]}/{{focal}}/{{IDs}}')
    benchmark:
        f'{config["path"]["root"]}/benchmarks/{{focal}}.{{IDs}}.crossMap.benchmark.txt'
    message:
        """
        This rule is an alternative implementation of the rule crossMap.
        Instead of taking each focal sample as a job and cross mapping in series using a for loop,
        here the cross mapping is done completely in parallel. 
        This implementation is not recommended, as it wastefully recreates a bwa index for each mapping
        operation. Use crossMap for smaller datasets or crossMap3 for larger datasets.
        """
    shell:
        """
        set +u;source activate {config[envs][metabagpipes]};set -u;
        cd $SCRATCHDIR

        echo -e "\nCopying assembly {input.contigs} and reads {input.R1} {input.R2} to $SCRATCHDIR"
        cp {input.contigs} {input.R1} {input.R2} .

        mkdir -p {output}

        # Define the focal sample ID, fsample: 
        # The one sample's assembly that reads will be mapped against
        fsampleID=$(echo $(basename $(dirname {input.contigs})))
        echo -e "\nFocal sample: $fsampleID ... "

        echo "Renaming and unzipping assembly ... "
        mv $(basename {input.contigs}) $(echo $fsampleID|sed 's/$/.fa.gz/g')
        gunzip $(echo $fsampleID|sed 's/$/.fa.gz/g')

        echo -e "\nIndexing assembly ... "
        bwa index $fsampleID.fa

        id=$(basename {output})
        echo -e "\nMapping reads from sample $id against assembly of focal sample $fsampleID ..."
        bwa mem -t {config[cores][crossMap]} $fsampleID.fa *.fastq.gz > $id.sam

        echo -e "\nDeleting no-longer-needed fastq files ... "
        rm *.gz
                
        echo -e "\nConverting SAM to BAM with samtools view ... " 
        samtools view -@ {config[cores][crossMap]} -Sb $id.sam > $id.bam

        echo -e "\nDeleting no-longer-needed sam file ... "
        rm $id.sam

        echo -e "\nSorting BAM file with samtools sort ... " 
        samtools sort -@ {config[cores][crossMap]} -o $id.sort $id.bam

        echo -e "\nDeleting no-longer-needed bam file ... "
        rm $id.bam

        echo -e "\nIndexing sorted BAM file with samtools index for CONCOCT input table generation ... " 
        samtools index $id.sort

        echo -e "\nCutting up assembly contigs >= 20kbp into 10kbp chunks and creating bedfile ... "
        cut_up_fasta.py $fsampleID.fa -c 10000 -o 0 --merge_last -b contigs_10K.bed > contigs_10K.fa

        echo -e "\nGenerating CONCOCT individual/intermediate coverage table ... "
        concoct_coverage_table.py contigs_10K.bed *.sort > ${{fsampleID}}_${{id}}_individual.tsv

        echo -e "\nCompressing CONCOCT coverage table ... "
        gzip ${{fsampleID}}_${{id}}_individual.tsv

        echo -e "\nRunning jgi_summarize_bam_contig_depths script to generate contig abundance/depth file for maxbin2 input ... "
        jgi_summarize_bam_contig_depths --outputDepth ${{fsampleID}}_${{id}}_individual.depth $id.sort

        echo -e "\nCompressing maxbin2/metabat2 depth file ... "
        gzip ${{fsampleID}}_${{id}}_individual.depth

        echo -e "\nMoving relevant files to {output}"
        mv *.gz {output}

        """

rule gatherCrossMap2: 
    input:
        expand(f'{config["path"]["root"]}/{config["folder"]["crossMap"]}/{{focal}}/{{IDs}}', focal = focal , IDs = IDs)
    shell:
        """
        echo idk
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

rule prepareRoary:
    input:
        taxonomy = rules.GTDBtkVis.output.text,
        binning = rules.binningVis.output.text,
        script = f'{config["path"]["root"]}/{config["folder"]["scripts"]}/{config["scripts"]["prepRoary"]}'
    output:
        directory(f'{config["path"]["root"]}/{config["folder"]["pangenome"]}/speciesBinIDs')
    benchmark:
        f'{config["path"]["root"]}/benchmarks/prepareRoary.benchmark.txt'
    message:
        """
        This rule matches the results from classifyGenomes->taxonomyVis with the completeness & contamination
        CheckM results from the metaWRAP reassembly->binningVis results, identifies speceies represented by 
        at least 10 high quality MAGs (completeness >= 90 & contamination <= 10), and outputs text files 
        with bin IDs for each such species. Also organizes the prokka output folders based on taxonomy.
        Note: Do not run this before finishing all prokka jobs!
        """
    shell:
        """
        set +u;source activate {config[envs][metabagpipes]};set -u
        cd $(dirname {input.taxonomy})

        echo -e "\nCreating speciesBinIDs folder containing.txt files with binIDs for each species that is represented by at least 10 high quality MAGs (completeness >= 90 & contamination <= 10) ... "
        Rscript {input.script}

        nSpecies=$(ls $(basename {output})|wc -l)
        nSpeciesTot=$(cat $(basename {output})/*|wc -l)
        nMAGsTot=$(paste {input.binning}|wc -l)
        echo -e "\nIdentified $nSpecies species represented by at least 10 high quality MAGs, totaling $nSpeciesTot MAGs out of $nMAGsTot total MAGs generated ... "

        echo -e "\nMoving speciesBinIDs folder to pangenome directory: $(dirname {output})"
        mv $(basename {output}) $(dirname {output})

        echo -e "\nOrganizing prokka folder according to taxonomy ... "
        echo -e "\nGFF files of identified species with at least 10 HQ MAGs will be copied to prokka/organzied/speciesSubfolder for roary input ... "
        cd $(dirname {output})
        mkdir -p prokka/organized

        for species in speciesBinIDs/*.txt;do

            speciesID=$(echo $(basename $species)|sed 's/.txt//g');
            echo -e "\nCreating folder and organizing prokka output for species $speciesID ... "
            mkdir -p prokka/organized/$speciesID

            while read line;do
                
                binID=$(echo $line|sed 's/.bin/_bin/g')
                echo "Copying GFF prokka output of bin $binID"
                cp prokka/unorganized/$binID/*.gff prokka/organized/$speciesID/

            done< $species
        done

        echo -e "\nDone"
        """ 


rule prepareRoaryMOTUS2:
    input:
        taxonomy = rules.taxonomyVis.output.text,
        binning = rules.binningVis.output.text,
        script = f'{config["path"]["root"]}/{config["folder"]["scripts"]}/{config["scripts"]["prepRoary"]}'
    output:
        directory(f'{config["path"]["root"]}/{config["folder"]["pangenome"]}/speciesBinIDs')
    benchmark:
        f'{config["path"]["root"]}/benchmarks/prepareRoary.benchmark.txt'
    message:
        """
        This rule matches the results from classifyGenomes->taxonomyVis with the completeness & contamination
        CheckM results from the metaWRAP reassembly->binningVis results, identifies speceies represented by 
        at least 10 high quality MAGs (completeness >= 90 & contamination <= 10), and outputs text files 
        with bin IDs for each such species. Also organizes the prokka output folders based on taxonomy.
        Note: Do not run this before finishing all prokka jobs!
        """
    shell:
        """
        set +u;source activate {config[envs][metabagpipes]};set -u
        cd $(dirname {input.taxonomy})

        echo -e "\nCreating speciesBinIDs folder containing.txt files with binIDs for each species that is represented by at least 10 high quality MAGs (completeness >= 90 & contamination <= 10) ... "
        Rscript {input.script}

        nSpecies=$(ls $(basename {output})|wc -l)
        nSpeciesTot=$(cat $(basename {output})/*|wc -l)
        nMAGsTot=$(paste {input.binning}|wc -l)
        echo -e "\nIdentified $nSpecies species represented by at least 10 high quality MAGs, totaling $nSpeciesTot MAGs out of $nMAGsTot total MAGs generated ... "

        echo -e "\nMoving speciesBinIDs folder to pangenome directory: $(dirname {output})"
        mv $(basename {output}) $(dirname {output})

        echo -e "\nOrganizing prokka folder according to taxonomy ... "
        echo -e "\nGFF files of identified species with at least 10 HQ MAGs will be copied to prokka/organzied/speciesSubfolder for roary input ... "
        cd $(dirname {output})
        mkdir -p prokka/organized

        for species in speciesBinIDs/*.txt;do

            speciesID=$(echo $(basename $species)|sed 's/.txt//g');
            echo -e "\nCreating folder and organizing prokka output for species $speciesID ... "
            mkdir -p prokka/organized/$speciesID

            while read line;do
                
                binID=$(echo $line|sed 's/.bin/_bin/g')
                echo "Copying GFF prokka output of bin $binID"
                cp prokka/unorganized/$binID/*.gff prokka/organized/$speciesID/

            done< $species
        done

        echo -e "\nDone"
        """ 

rule roaryTop10:
    input:
        f'{config["path"]["root"]}/{config["folder"]["pangenome"]}/prokka/organized/'
    output:
        directory(f'{config["path"]["root"]}/{config["folder"]["pangenome"]}/roary/top10/')
    benchmark:
        f'{config["path"]["root"]}/benchmarks/roaryTop10.roary.benchmark.txt'
    message:
        """
        Runs pangenome for ~692 MAGs belonging to 10 species:
        Agathobacter rectale, Bacteroides uniformis, 
        Ruminococcus_E bromii_B, Gemmiger sp003476825, 
        Blautia_A wexlerae, Dialister invisus,
        Anaerostipes hadrus, Fusicatenibacter saccharivorans,
        Eubacterium_E hallii, and NA
        """
    shell:
        """
        set +u;source activate prokkaroary;set -u
        mkdir -p $(dirname {output})
        cd $SCRATCHDIR

        cp -r {input}/Agathobacter_rectale/* . 
        cp -r {input}/Bacteroides_uniformis/* . 
        cp -r {input}/Ruminococcus_E_bromii_B/* . 
        cp -r {input}/Gemmiger_sp003476825/* . 
        cp -r {input}/Blautia_A_wexlerae/* .
        cp -r {input}/Dialister_invisus/* . 
        cp -r {input}/Anaerostipes_hadrus/* . 
        cp -r {input}/Fusicatenibacter_saccharivorans/* . 
        cp -r {input}/Eubacterium_E_hallii/* . 
        cp -r {input}/NA/* .
                
        roary -s -p {config[cores][roary]} -i {config[params][roaryI]} -cd {config[params][roaryCD]} -f yes_al -e -v *.gff
        cd yes_al
        create_pan_genome_plots.R 
        cd ..
        mkdir -p {output}

        mv yes_al/* {output}
        """

rule phylophlan:
    input:
        f'/home/zorrilla/workspace/european/dna_bins'
    output:
        directory(f'/scratch/zorrilla/phlan/out')
    benchmark:
        f'/scratch/zorrilla/phlan/logs/bench.txt'
    shell:
        """
        cd $SCRATCHDIR
        cp -r {input} . 
        cp $(dirname {output})/*.cfg .
        mkdir -p logs

        phylophlan -i dna_bins \
                    -d phylophlan \
                    -f 02_tol.cfg \
                    --genome_extension fa \
                    --diversity low \
                    --fast \
                    -o out \
                    --nproc 128 \
                    --verbose 2>&1 | tee logs/phylophlan.logs

        cp -r out $(dirname {output})
        """

rule phylophlanPlant:
    input:
        f'/home/zorrilla/workspace/china_soil/dna_bins'
    output:
        directory(f'/home/zorrilla/workspace/china_soil/phlan/')
    benchmark:
        f'/scratch/zorrilla/phlan/logs/benchPlant.txt'
    shell:
        """
        cd $SCRATCHDIR
        cp -r {input} . 
        cp /scratch/zorrilla/phlan/*.cfg .
        mkdir -p logs

        phylophlan -i dna_bins \
                    -d phylophlan \
                    -f 02_tol.cfg \
                    --genome_extension fa \
                    --diversity low \
                    --fast \
                    -o $(basename {output}) \
                    --nproc 128 \
                    --verbose 2>&1 | tee logs/phylophlan.logs

        cp -r $(basename {output}) $(dirname {output})
        """
        
rule phylophlanMeta:
    input:
        f'/home/zorrilla/workspace/european/dna_bins'
    output:
        directory(f'/home/zorrilla/workspace/european/phlan/dist')
    benchmark:
        f'/scratch/zorrilla/phlan/logs/bench.txt'
    shell:
        """
        cd {input}
        cd ../

        phylophlan_metagenomic -i $(basename {input}) -o $(basename {output})_dist --nproc 2 --only_input
    
        mv $(basename {output})_dist $(basename {output})
        mv -r $(basename {output}) $(dirname {output})
        """


rule phylophlanMetaAll:
    input:
        lab=f'/home/zorrilla/workspace/korem/dna_bins',
        gut=f'/home/zorrilla/workspace/european/dna_bins' ,
        plant=f'/home/zorrilla/workspace/china_soil/dna_bins' ,
        soil=f'/home/zorrilla/workspace/straya/dna_bins' ,
        ocean=f'/scratch/zorrilla/dna_bins' 
    output:
        directory(f'/home/zorrilla/workspace/european/phlan/all')
    benchmark:
        f'/scratch/zorrilla/phlan/logs/allMetaBench.txt'
    shell:
        """
        mkdir -p {output}
        cd $SCRATCHDIR

        mkdir allMAGs
        cp {input.lab}/* allMAGs
        cp {input.gut}/* allMAGs
        cp {input.plant}/* allMAGs
        cp {input.soil}/* allMAGs
        cp {input.ocean}/* allMAGs

        phylophlan_metagenomic -i allMAGs -o all --nproc 4 --only_input
        mv all_distmat.tsv $(dirname {output})
        """

rule drawTree:
    input:
        f'/home/zorrilla/workspace/china_soil/phlan'
    shell:
        """
        cd {input}
        graphlan.py dna_bins.tre.iqtree tree.out
        """

rule makePCA:
    input:
        f'/home/zorrilla/workspace/european/phlan'
    shell:
        """
        cd $SCRATCHDIR
        echo -e "\nCopying files to scratch dir: $SCRATCHDIR"
        cp {input}/*.tsv {input}/*.ids {input}/*.R .

        echo -e "\nRunning nmds.R script ... "
        Rscript nmds.R
        rm *.tsv *.ids *.R

        mkdir -p nmds
        mv *.pdf nmds
        mv nmds {input}
        """

rule drep:
    input:
        f'{config["path"]["root"]}/dna_bins'
    output:
        directory(f'{config["path"]["root"]}/drep_drep')
    benchmark:
        f'{config["path"]["root"]}/benchmarks/drep_drep.benchmark.txt'
    shell:
        """
        set +u;source activate drep;set -u
        cp -r {input} $SCRATCHDIR
        cd $SCRATCHDIR

        dRep dereplicate drep_drep -g $(basename {input})/*.fa -p 48 -comp 50 -con 10
        mv drep_drep $(dirname {input})
        """

rule drepComp:
    input:
        f'{config["path"]["root"]}/dna_bins'
    output:
        directory(f'{config["path"]["root"]}/drep_comp')
    benchmark:
        f'{config["path"]["root"]}/benchmarks/drep_comp.benchmark.txt'
    shell:
        """
        set +u;source activate drep;set -u
        cp -r {input} $SCRATCHDIR
        cd $SCRATCHDIR

        dRep compare drep_comp -g $(basename {input})/*.fa -p 48
        mv drep_comp $(dirname {input})
        """
