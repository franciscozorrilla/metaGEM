rule kallisto2concoctTable: 
    input:
        f'{config["path"]["root"]}/{config["folder"]["kallisto"]}/{{focal}}/'
    output: 
        f'{config["path"]["root"]}/{config["folder"]["concoct"]}/{{focal}}/cov/coverage_table.tsv'
    message:
        """
        This rule is necessary for the crossMapParallel implementation subworkflow.
        It summarizes the individual concoct input tables for a given focal sample.
        Note: silence output if not using parallel mapping approach
        """
    shell:
        """
        # Activate metagem environment
        set +u;source activate {config[envs][metagem]};set -u;

        # Create output folder
        mkdir -p $(dirname {output})

        # Compile individual mapping results into coverage table for given assembly
        python {config[path][root]}/{config[folder][scripts]}/{config[scripts][kallisto2concoct]} \
            --samplenames <(for s in {input}/*; do echo $s|sed 's|^.*/||'; done) \
            $(find {input} -name "*.gz") > {output}
    
        """