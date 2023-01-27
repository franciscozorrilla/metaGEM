rule metabat_single:
    input:
        assembly = rules.megahit.output,
        R1 = rules.qfilter.output.R1,
        R2 = rules.qfilter.output.R2
    output:
        directory(f'{config["path"]["root"]}/{config["folder"]["metabat"]}/{{IDs}}/{{IDs}}.metabat-bins')
    benchmark:
        f'{config["path"]["root"]}/{config["folder"]["benchmarks"]}/{{IDs}}.metabat.benchmark.txt'
    message:
        """
        Implementation of metabat2 where only coverage information from the focal sample is used
        for binning. Use with the crossMapParallel subworkflow, where cross sample coverage information
        is only used by CONCOCT.
        """
    shell:
        """
        # Activate metagem environment
        set +u;source activate {config[envs][metagem]};set -u;
        # Make job specific scratch dir
        fsampleID=$(echo $(basename $(dirname {input.assembly})))
        echo -e "\nCreating temporary directory {config[path][scratch]}/{config[folder][metabat]}/${{fsampleID}} ... "
        mkdir -p {config[path][scratch]}/{config[folder][metabat]}/${{fsampleID}}

        # Move into scratch dir
        cd {config[path][scratch]}/{config[folder][metabat]}/${{fsampleID}}

        # Copy files
        cp {input.assembly} {input.R1} {input.R2} .

        echo -e "\nFocal sample: $fsampleID ... "

        echo "Renaming and unzipping assembly ... "
        mv $(basename {input.assembly}) $(echo $fsampleID|sed 's/$/.fa.gz/g')
        gunzip $(echo $fsampleID|sed 's/$/.fa.gz/g')

        echo -e "\nIndexing assembly ... "
        bwa index $fsampleID.fa

        id=$(basename {output})
        echo -e "\nMapping reads from sample against assembly $fsampleID ..."
        bwa mem -t {config[cores][metabat]} $fsampleID.fa *.fastq.gz > $id.sam

        echo -e "\nDeleting no-longer-needed fastq files ... "
        rm *.gz
                
        echo -e "\nConverting SAM to BAM with samtools view ... " 
        samtools view -@ {config[cores][metabat]} -Sb $id.sam > $id.bam

        echo -e "\nDeleting no-longer-needed sam file ... "
        rm $id.sam

        echo -e "\nSorting BAM file with samtools sort ... " 
        samtools sort -@ {config[cores][metabat]} -o $id.sort $id.bam

        echo -e "\nDeleting no-longer-needed bam file ... "
        rm $id.bam

        # Run metabat2
        echo -e "\nRunning metabat2 ... "
        jgi_summarize_bam_contig_depths --outputDepth $id.depth.txt $id.sort

        metabat2 -i $fsampleID.fa -a $id.depth.txt -s \
            {config[params][metabatMin]} \
            -v --seed {config[params][seed]} \
            -t 0 -m {config[params][minBin]} \
            -o $(basename $(dirname {output}))

        rm $fsampleID.fa
        rm $id.depth.txt

        # Move files to output dir
        mv *.fa {output}
        """