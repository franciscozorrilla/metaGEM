rule maxbin_single:
    input:
        assembly = rules.megahit.output,
        R1 = rules.qfilter.output.R1,
        R2 = rules.qfilter.output.R2
    output:
        directory(f'{config["path"]["root"]}/{config["folder"]["maxbin"]}/{{IDs}}/{{IDs}}.maxbin-bins')
    benchmark:
        f'{config["path"]["root"]}/{config["folder"]["benchmarks"]}/{{IDs}}.maxbin.benchmark.txt'
    shell:
        """
        # Activate metagem environment
        set +u;source activate {config[envs][metagem]};set -u;

        # Create output folder
        mkdir -p $(dirname {output})

        # Make job specific scratch dir
        fsampleID=$(echo $(basename $(dirname {input.assembly})))
        echo -e "\nCreating temporary directory {config[path][scratch]}/{config[folder][maxbin]}/${{fsampleID}} ... "
        mkdir -p {config[path][scratch]}/{config[folder][maxbin]}/${{fsampleID}}

        # Move into scratch dir
        cd {config[path][scratch]}/{config[folder][maxbin]}/${{fsampleID}}

        # Copy files to tmp
        cp -r {input.assembly} {input.R1} {input.R2} .

        echo -e "\nUnzipping assembly ... "
        gunzip $(basename {input.assembly})
        
        echo -e "\nRunning maxbin2 ... "
        run_MaxBin.pl -contig $(echo $(basename {input.assembly})|sed 's/.gz//') \
            -out $(basename $(dirname {output})) \
            -reads $(basename {input.R1}) \
            -reads2 $(basename {input.R2}) \
            -thread {config[cores][maxbin]}

        rm $(echo $(basename {input.assembly})|sed 's/.gz//')
        
        mkdir $(basename {output})
        mv *.fasta *.summary *.abundance *.abund1 *.abund2 $(basename {output})
        mv $(basename {output}) $(dirname {output})

        """