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
