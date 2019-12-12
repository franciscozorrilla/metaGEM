configfile: "config.yaml"

import os
import glob

def get_ids_from_path_pattern(path_pattern):
    ids = sorted([os.path.basename(os.path.splitext(val)[0])
                  for val in (glob.glob(path_pattern))])
    return ids

# Make sure that final_bins/ folder contains all bins in single folder for binIDs
# wildcard to work. Use moveBins rule or perform manually.
binIDs = get_ids_from_path_pattern('final_bins/*.faa')
IDs = get_ids_from_path_pattern('dataset/*')

DATA_READS_1 = f'{config["path"]["root"]}/{config["folder"]["data"]}/{{IDs}}/{{IDs}}_1.fastq.gz'
DATA_READS_2 = f'{config["path"]["root"]}/{config["folder"]["data"]}/{{IDs}}/{{IDs}}_2.fastq.gz'


rule all:
    input:
        expand(f'{config["path"]["root"]}/{config["folder"]["assemblies"]}/{{IDs}}/contigs.fasta.gz', IDs=IDs)
    shell:
        """
        echo {input}
        """


rule createFolders:
    input:
        config["path"]["root"]
    message:
        """
        Not strictly necessary, as we could create individual folders within rules, 
        but doing this at the beginning may give the user a better idea of what 
        output/intermediate files to expect
        """
    shell:
        """
        cd {input}
        paste config.yaml |cut -d':' -f2|tail -n +4|head -n 17 > folders.txt
        while read line;do mkdir -p $line;done < folders.txt
        rm folders.txt
        """


rule downloadToy:
    input:
        config["dbs"]["toy"]
    shell:
        """
        cd {config[path][root]}/{config[folder][data]}
        while read line;do wget $line;done < {input}
        for file in *;do mv $file ./$(echo $file|sed 's/?download=1//g');done
        """


rule organizeData:
    input:
        config["path"]["root"]
    message:
        """
        Assuming all samples are downloaded to the same directory, 
        sorting paired end raw reads into sample specific sub folders within the 'dataset' folder.
        Note: Use with caution, may malfunction in corner cases.
        """
    shell:
        """
        cd {input}/{config[folder][data]}
        for file in *.gz; do echo $file; done | sed 's/_.*$//g' | sed 's/.fastq.gz//g' | uniq > ID_samples.txt
        while read line; do mkdir -p $line; mv $line*.gz $line; done < ID_samples.txt
        rm ID_samples.txt
        """


rule qfilter: 
    input:
        R1 = DATA_READS_1,
        R2 = DATA_READS_2
    output:
        R1 = config["path"]["root"]+"/"+config["folder"]["qfiltered"]+"/{IDs}/{IDs}_1.fastq.gz", 
        R2 = config["path"]["root"]+"/"+config["folder"]["qfiltered"]+"/{IDs}/{IDs}_2.fastq.gz" 
    shell:
        """
        set +u;source activate {config[envs][metabagpipes]};set -u;
        mkdir -p $(dirname $(dirname {output.R1}))
        mkdir -p $(dirname {output.R1})
        fastp --thread {config[cores][fastp]} -i {input.R1} -I {input.R2} -o {output.R1} -O {output.R2} -j $(dirname {output.R1})/$(echo $(basename $(dirname {output.R1}))).json -h $(dirname {output.R1})/$(echo $(basename $(dirname {output.R1}))).html
        """

rule megahit:
    input:
        R1 = rules.qfilter.output.R1, 
        R2 = rules.qfilter.output.R2
    output:
        config["path"]["root"]+"/"+config["folder"]["assemblies"]+"/{IDs}/contigs.fasta.gz"
    benchmark:
        config["path"]["root"]+"/"+"benchmarks/{IDs}.megahit.benchmark.txt"
    shell:
        """
        set +u;source activate {config[envs][metabagpipes]};set -u;
        mkdir -p $(dirname {output})
        cd $TMPDIR
        cp {input.R1} {input.R2} $TMPDIR
        megahit -t {config[cores][megahit]} --presets meta-large --verbose -1 $(basename {input.R1}) -2 $(basename {input.R2}) -o tmp
        mv tmp/final.contigs.fa contigs.fasta
        gzip contigs.fasta
        mv contigs.fasta.gz $(dirname {output})
        """


rule assemblyVisMetaspades:
    input:
        config["path"]["root"]
    message:
        """
        This may take a long time to run with many samples. Run the lines below to run rule in background of login node instead of through snakemake.
        Raw data: nohup sh -c 'for folder in */;do echo -n "$folder "|sed "s|/||g" >> dataset.stats;zcat "$folder"*_1.fastq.gz | awk "{{if(NR%4==2) print length($1)}}" | sort | uniq -c >> dataset.stats;done' &
        Assembly data: nohup sh -c 'for folder in */;do for file in $folder*.gz;do N=$(less $file|grep -c ">"); L=$(less $file|grep ">"|cut -d "_" -f4|awk '"'"'{sum+=$NF} END{print sum}'"'"'); C=$(less $file|grep ">"|cut -d "_" -f6|awk '"'"'{sum+=$NF} END { if (NR > 0) print sum / NR }'"'"'); echo $(echo $file|sed "s|/contigs.fasta.gz||g") $N $L $C >> assembly.stats; done;done'&
        """
    shell:
        """
        set +u;source activate memotenv;set -u;
        cd {input}/{config[folder][data]}
        for folder in */;do 
        echo -n "$folder "|sed "s|/||g" >> dataset.stats;
        zcat "$folder"*_1.fastq.gz | awk '{{if(NR%4==2) print length($1)}}' | sort | uniq -c >> dataset.stats;
        done
        mv dataset.stats {input}/{config[folder][assemblies]}
        cd {input}/{config[folder][assemblies]}
        for folder in */;do 
        for file in $folder*.gz;do 
        N=$(less $file|grep -c ">"); 
        L=$(less $file|grep ">"|cut -d '_' -f4|awk '{{sum+=$1}} END{{print sum}}');
        A=$(awk -v n="$N" -v l="$L" 'BEGIN{{ if (n>0) print l / n}}') 
        C=$(less $file|grep ">"|cut -d '_' -f6|awk '{{sum+=$1}} END {{ if (NR > 0) print sum / NR }}'); 
        echo $(echo $file|sed 's|/contigs.fasta.gz||g') $N $L $A $C >> assembly.stats;
        done;
        done
        Rscript {config[path][root]}/{config[folder][scripts]}/{config[scripts][assemblyVis]}
        """

rule assemblyVisMegahit:
    shell:
        """
        for folder in */;do 
            for file in $folder*.gz;do 
                N=$(less $file|grep -c ">"); 
                L=$(less $file|grep ">"|cut -d ' ' -f4|sed 's/len=//'|awk '{{sum+=$1}}END{{print sum}}');
                A=$(awk -v n="$N" -v l="$L" 'BEGIN{{ if (n>0) print l / n}}'); 
                M=$(less $file|grep ">"|cut -d ' ' -f4|sed 's/len=//'|sort -n | awk 'NF{{a[NR]=$1;c++}}END{{print (c%2==0)?((a[c/2]+a[c/2+1])/2):a[c/2+1]}}');
                T=$(less $file|grep ">"|cut -d ' ' -f4|sed 's/len=//'|awk '$1>=1000{c++} END{print c+0}');
                S=$(less $file|grep ">"|cut -d ' ' -f4|sed 's/len=//'|awk '$1>=1000'|awk '{{sum+=$1}}END{{print sum}}')
                echo $(echo $file|sed 's|/contigs.fasta.gz||g') $N $L $A $M $T $S>> assembly.stats;
            done;
        done
        """

rule kallisto:
    input:
        contigs = rules.metaspades.output,
        reads = f'{config["path"]["root"]}/{config["folder"]["qfiltered"]}/'
    output:
        f'{config["path"]["root"]}/{config["folder"]["concoctInput"]}/{{IDs}}_concoct_inputtableR.tsv'
    benchmark:
        f'{config["path"]["root"]}/benchmarks/{{IDs}}.kallisto.benchmark.txt'
    shell:
        """
        set +u;source activate {config[envs][metabagpipes]};set -u;
        cp {input.contigs} $TMPDIR
        cd $TMPDIR
        gunzip $(basename {input.contigs})
        cut_up_fasta.py -c {config[params][cutfasta]} -o 0 -m contigs.fasta > metaspades_c10k.fa
        
        kallisto index metaspades_c10k.fa -i $(basename $(dirname {input.contigs})).kaix
        
        for folder in {input.reads}*;do
        cp $folder/*.fastq.gz $TMPDIR;
        
        kallisto quant --threads {config[cores][kallisto]} --plaintext \
            -i $(basename $(dirname {input.contigs})).kaix \
            -o . *_1.fastq.gz *_2.fastq.gz;
        
        gzip abundance.tsv;
        mv abundance.tsv.gz $(echo $(basename $folder)_abundance.tsv.gz);
        rm *.fastq.gz;
        done
        
        python {config[scripts][kallisto2concoct]} \
            --samplenames <(for s in *abundance.tsv.gz; do echo $s | sed 's/_abundance.tsv.gz//'g; done) *abundance.tsv.gz > $(basename {output})
        
        mv $(basename {output}) $(dirname {output})
        """


rule concoct:
    input:
        table = rules.kallisto.output,
        contigs = rules.metaspades.output
    output:
        directory(f'{config["path"]["root"]}/{config["folder"]["concoctOutput"]}/{{IDs}}/{{IDs}}.concoct-bins')
    benchmark:
        f'{config["path"]["root"]}/benchmarks/{{IDs}}.concoct.benchmark.txt'
    shell:
        """
        set +u;source activate {config[envs][metabagpipes]};set -u;
        mkdir -p $(dirname {output})
        mkdir -p $(dirname $(dirname {output}))
        cp {input.contigs} {input.table} $TMPDIR
        cd $TMPDIR
        gunzip $(basename {input.contigs})
        cut_up_fasta.py -c {config[params][cutfasta]} -o 0 -m contigs.fasta > metaspades_c10k.fa
        
        concoct --coverage_file {input.table} --composition_file metaspades_c10k.fa \
            -b $(basename $(dirname {output})) \
            -t {config[cores][concoct]} \
            -c {config[params][concoct]}
            
        merge_cutup_clustering.py $(basename $(dirname {output}))_clustering_gt1000.csv > $(basename $(dirname {output}))_clustering_merged.csv
        
        mkdir -p $(basename {output})
        extract_fasta_bins.py contigs.fasta $(basename $(dirname {output}))_clustering_merged.csv --output_path $(basename {output})
        mv $(basename {output}) *.log *.txt *.csv *.tab $(dirname {output})
        """


rule metabat:
    input:
        assembly = rules.metaspades.output,
        R1 = rules.qfilter.output.R1, 
        R2 = rules.qfilter.output.R2
    output:
        directory(f'{config["path"]["root"]}/{config["folder"]["metabat"]}/{{IDs}}/{{IDs}}.metabat-bins')
    benchmark:
        f'{config["path"]["root"]}/benchmarks/{{IDs}}.metabat.benchmark.txt'
    shell:
        """
        set +u;source activate {config[envs][metabagpipes]};set -u;
        mkdir -p $(dirname $(dirname {output}))
        mkdir -p $(dirname {output})
        cp {input.assembly} {input.R1} {input.R2} $TMPDIR
        cd $TMPDIR
        mv $(basename {input.assembly}) $(basename $(dirname {input.assembly})).gz
        gunzip $(basename $(dirname {input.assembly})).gz
        bwa index $(basename $(dirname {input.assembly}))
        
        bwa mem -t {config[cores][metabat]} $(basename $(dirname {input.assembly})) \
            $(basename {input.R1}) \
            $(basename {input.R2}) > $(basename $(dirname {input.assembly})).sam
         
        samtools view -@ {config[cores][metabat]} -Sb $(basename $(dirname {input.assembly})).sam > $(basename $(dirname {input.assembly})).bam
        samtools sort -@ {config[cores][metabat]} $(basename $(dirname {input.assembly})).bam > $(basename $(dirname {input.assembly})).sort
        runMetaBat.sh $(basename $(dirname {input.assembly})) $(basename $(dirname {input.assembly})).sort
        mv *.txt *.tab $(basename {output}) $(dirname {output})
        """

rule maxbin:
    input:
        assembly = rules.metaspades.output,
        R1 = rules.qfilter.output.R1, 
        R2 = rules.qfilter.output.R2
    output:
        directory(f'{config["path"]["root"]}/{config["folder"]["maxbin"]}/{{IDs}}/{{IDs}}.maxbin-bins')
    benchmark:
        f'{config["path"]["root"]}/benchmarks/{{IDs}}.maxbin.benchmark.txt'
    shell:
        """
        set +u;source activate {config[envs][metabagpipes]};set -u;
        mkdir -p $(dirname $(dirname {output}))
        mkdir -p $(dirname {output})
        cp {input.assembly} {input.R1} {input.R2} $TMPDIR
        cd $TMPDIR
        gunzip contigs.fasta.gz
        
        run_MaxBin.pl -contig contigs.fasta -out $(basename $(dirname {output})) \
            -reads *_1.fastq.gz -reads2 *_2.fastq.gz \
            -thread {config[cores][maxbin]}
        
        rm contigs.fasta *.gz
        mkdir $(basename {output})
        mv *.fasta $(basename {output})
        mv *.abund1 *.abund2 *.tab $(basename {output}) $(dirname {output})
        """


rule binRefine:
    input:
        concoct = f'{config["path"]["root"]}/{config["folder"]["concoctOutput"]}/{{IDs}}/{{IDs}}.concoct-bins',
        metabat = f'{config["path"]["root"]}/{config["folder"]["metabat"]}/{{IDs}}/{{IDs}}.metabat-bins',
        maxbin = f'{config["path"]["root"]}/{config["folder"]["maxbin"]}/{{IDs}}/{{IDs}}.maxbin-bins'
    output:
        directory(f'{config["path"]["root"]}/{config["folder"]["refined"]}/{{IDs}}')
    benchmark:
        f'{config["path"]["root"]}/benchmarks/{{IDs}}.binRefine.benchmark.txt'
    shell:
        """
        set +u;source activate {config[envs][metawrap]};set -u;
        cp -r {input.concoct} {input.metabat} {input.maxbin} $TMPDIR
        mkdir -p $(dirname {output})
        cd $TMPDIR
        
        metaWRAP bin_refinement -o . \
            -A $(basename {input.concoct}) \
            -B $(basename {input.metabat}) \
            -C $(basename {input.maxbin}) \
            -t {config[cores][refine]} \
            -m {config[params][refineMem]} \
            -c {config[params][refineComp]} \
            -x {config[params][refineCont]}
 
        rm -r $(basename {input.concoct}) $(basename {input.metabat}) $(basename {input.maxbin})
        mkdir -p {output}
        mv * {output}
        """


rule binReassemble:
    input:
        R1 = rules.qfilter.output.R1, 
        R2 = rules.qfilter.output.R2,
        refinedBins = f'{config["path"]["root"]}/{config["folder"]["refined"]}/{{IDs}}/metawrap_50_10_bins'
    output:
        directory(f'{config["path"]["root"]}/{config["folder"]["reassembled"]}/{{IDs}}')
    benchmark:
        f'{config["path"]["root"]}/benchmarks/{{IDs}}.binReassemble.benchmark.txt'
    shell:
        """
        set +u;source activate {config[envs][metawrap]};set -u;
        mkdir -p $(dirname {output})
        cp -r {input.refinedBins} {input.R1} {input.R2} $TMPDIR
        cd $TMPDIR
        
        metaWRAP reassemble_bins -o $(basename {output}) \
            -b $(basename {input.refinedBins}) \
            -1 $(basename {input.R1}) \
            -2 $(basename {input.R1}) \
            -t {config[cores][reassemble]} \
            -m {config[params][reassembleMem]} \
            -c {config[params][reassembleComp]} \
            -x {config[params][reassembleCont]}
        
        rm -r $(basename {input.refinedBins})
        rm -r $(basename {output})/work_files
        rm *.fastq.gz 
        mv * $(dirname {output})
        """


rule binningVis:
    input:
        config["path"]["root"]
    message:
        """
        Generate bar plot with number of bins and density plot of bin contigs, 
        total length, completeness, and contamination across different tools.
        """
    shell:
        """
        set +u;source activate memotenv;set -u;
        cd {input}/{config[folder][concoctOutput]}
        for folder in */;do 
            var=$(echo $folder|sed 's|/||g'); 
            for bin in $folder*concoct-bins/*.fa;do 
                name=$(echo $bin | sed "s|^.*/|$var.bin.|g" | sed 's/.fa//g'); 
                N=$(less $bin | grep -c ">");
                C=$(less $bin | grep ">" | cut -d '_' -f6 | awk '{{sum+=$1}} END {{ if (NR > 0) print sum / NR }}');
                echo $name $N $C >> concoct_bins.stats;
            done;
        done
        mv *.stats {input}/{config[folder][reassembled]}
        cd {input}/{config[folder][metabat]}
        for folder in */;do 
            var=$(echo $folder | sed 's|/||'); 
            for bin in $folder*metabat-bins/*.fa;do 
                name=$(echo $bin|sed 's/.fa//g'|sed 's|^.*/||g'|sed "s/^/$var./g"); 
                N=$(less $bin | grep -c ">");
                C=$(less $bin | grep ">" | cut -d '_' -f6 | awk '{{sum+=$1}} END {{ if (NR > 0) print sum / NR }}');
                echo $name $N $C >> metabat_bins.stats;
            done;
        done
        mv *.stats {input}/{config[folder][reassembled]}
        cd {input}/{config[folder][maxbin]}
        for folder in */;do
            for bin in $folder*maxbin-bins/*.fasta;do 
                name=$(echo $bin | sed 's/.fasta//g' | sed 's|^.*/||g');
                N=$(less $bin | grep -c ">");
                C=$(less $bin | grep ">"|cut -d '_' -f6 | awk '{{sum+=$1}} END {{ if (NR > 0) print sum / NR }}');
                echo $name $N $C >> maxbin_bins.stats;
            done;
        done
        mv *.stats {input}/{config[folder][reassembled]}
        cd {input}/{config[folder][refined]}
        for folder in */;do
            var=$(echo $folder | sed 's|/||g');paste $folder*concoct-bins.stats|tail -n +2 | sed "s/^/$var.bin./g";
        done >> concoct.checkm
        for folder in */;do 
            var=$(echo $folder | sed 's|/||g');paste $folder*metabat-bins.stats|tail -n +2 | sed "s/^/$var./g";
        done >> metabat.checkm
        for folder in */;do paste $folder*maxbin-bins.stats|tail -n +2;done >> maxbin.checkm
        for folder in */;do var=$(echo $folder|sed 's|/||g');paste $folder*etawrap_50_10_bins.stats|tail -n +2|sed "s/^/$var./g";done >> refined.checkm
        for folder in */;do 
            samp=$(echo $folder | sed 's|/||');
            for bin in $folder*metawrap_50_10_bins/*.fa;do 
                name=$(echo $bin | sed 's/.fa//g'|sed 's|^.*/||g'|sed "s/^/$samp./g");
                N=$(less $bin | grep -c ">");
                C=$(less $bin | grep ">"|cut -d '_' -f6|awk '{{sum+=$1}} END {{ if (NR > 0) print sum / NR }}');
                echo $name $N $C >> refined_bins.stats;
            done;
        done
        mv *.stats *.checkm {input}/{config[folder][reassembled]}
        cd {input}/{config[folder][reassembled]}
        for folder in */;do 
            samp=$(echo $folder | sed 's|/||');
            for bin in $folder*reassembled_bins/*.fa;do 
                name=$(echo $bin | sed 's/.fa//g' | sed 's|^.*/||g' | sed "s/^/$samp./g");
                N=$(less $bin | grep -c ">");
                C=$(less $bin | grep ">"|cut -d '_' -f6|awk '{{sum+=$1}} END {{ if (NR > 0) print sum / NR }}');
                echo $name $N $C >> reassembled_bins.stats;
            done;
        done
        for folder in */;do var=$(echo $folder|sed 's|/||g');paste $folder*reassembled_bins.stats|tail -n +2|sed "s/^/$var./g";done >> reassembled.checkm
        Rscript {config[path][root]}/{config[folder][scripts]}/{config[scripts][binningVis]}
        rm Rplots.pdf
        """


rule classifyGenomes:
    input:
        bins = f'{config["path"]["root"]}/{config["folder"]["reassembled"]}/{{IDs}}/reassembled_bins',
        script = f'{config["path"]["root"]}/{config["folder"]["scripts"]}/classify-genomes'
    output:
        directory(f'{config["path"]["root"]}/{config["folder"]["classification"]}/{{IDs}}')
    benchmark:
        f'{config["path"]["root"]}/benchmarks/{{IDs}}.classify-genomes.benchmark.txt'
    shell:
        """
        set +u;source activate {config[envs][metabagpipes]};set -u;
        mkdir -p {output}
        cd $TMPDIR
        cp -r {input.script}/* {input.bins}/* .
        for bin in *.fa; do
        echo "RUNNING BIN $bin"
        $PWD/classify-genomes $bin -t {config[cores][classify]} -o $(echo $bin|sed 's/.fa/.taxonomy/')
        echo "DONE" 
        echo " "
        cp *.taxonomy {output}
        rm *.taxonomy
        rm $bin 
        done
        """


rule abundance:
    input:
        bins = f'{config["path"]["root"]}/{config["folder"]["reassembled"]}/{{IDs}}/reassembled_bins',
        R1 = rules.qfilter.output.R1, 
        R2 = rules.qfilter.output.R2
    output:
        directory(f'{config["path"]["root"]}/{config["folder"]["abundance"]}/{{IDs}}')
    benchmark:
        f'{config["path"]["root"]}/benchmarks/{{IDs}}.abundance.benchmark.txt'
    message:
        """
        Calculate bin abundance fraction using the following:

        binAbundanceFraction = (L * X / Y / Z) * 100

        X = # of reads mapped to bin_i from sample_k
        Y = length of bin_i
        Z = # of reads mapped to all bins in sample_k
        L = length of reads (100 bp)
        """
    shell:
        """
        set +u;source activate {config[envs][metabagpipes]};set -u;
        mkdir -p {output}
        cd $TMPDIR
        cp {input.bins}/* .
        cat *.fa > $(basename {output}).fa
        cp {input.R1} {input.R2} .
        echo "CREATING INDEX FOR BIN CONCATENATION AND MAPPING READS "
        bwa index $(basename {output}).fa
        bwa mem -t {config[cores][abundance]} $(basename {output}).fa \
            $(basename {input.R1}) $(basename {input.R2}) > $(basename {output}).sam
        samtools view -@ {config[cores][abundance]} -Sb $(basename {output}).sam > $(basename {output}).bam
        samtools sort -@ {config[cores][abundance]} $(basename {output}).bam $(basename {output}).sort
        samtools flagstat $(basename {output}).sort.bam > map.stats
        cp map.stats {output}/$(basename {output})_map.stats
        rm $(basename {output}).fa
        echo "DONE MAPPING READS TO BIN CONCATENATION, BEGIN MAPPING READS TO EACH BIN "
        for bin in *.fa;do
            mkdir -p $(echo "$bin"| sed "s/.fa//")
            mv $bin $(echo "$bin"| sed "s/.fa//")
            cd $(echo "$bin"| sed "s/.fa//")
            echo "INDEXING AND MAPPING BIN $bin "
            bwa index $bin
            bwa mem -t {config[cores][abundance]} $bin \
                ../$(basename {input.R1}) ../$(basename {input.R2}) > $(echo "$bin"|sed "s/.fa/.sam/")
            samtools view -@ {config[cores][abundance]} -Sb $(echo "$bin"|sed "s/.fa/.sam/") > $(echo "$bin"|sed "s/.fa/.bam/")
            samtools sort -@ {config[cores][abundance]} $(echo "$bin"|sed "s/.fa/.bam/") $(echo "$bin"|sed "s/.fa/.sort/")
            samtools flagstat $(echo "$bin"|sed "s/.fa/.sort.bam/") > $(echo "$bin"|sed "s/.fa/.map/")
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
            rm -r $(echo "$bin"| sed "s/.fa//")
        done
        cat *.abund > $(basename {output}).abund
        mv $(basename {output}).abund {output}
        """


rule taxonomyVis:
    shell:
        """
        cd {config[path][root]}/{config[folder][classification]}
        for folder in */;do 
            for file in $folder*.taxonomy;do 
                fasta=$(echo $file | sed 's|^.*/||' | sed 's/.taxonomy//g' | sed 's/.orig//g' | sed 's/.permissive//g' | sed 's/.strict//g'); 
                NCBI=$(less $file | grep NCBI | cut -d ' ' -f4);
                tax=$(less $file | grep tax | sed 's/Consensus taxonomy: //g');
                motu=$(less $file | grep mOTUs | sed 's/Consensus mOTUs: //g');
                detect=$(less $file | grep detected | sed 's/Number of detected genes: //g');
                percent=$(less $file | grep agreeing | sed 's/Percentage of agreeing genes: //g' | sed 's/%//g');
                map=$(less $file | grep mapped | sed 's/Number of mapped genes: //g');
                cog=$(less $file | grep COG | cut -d$'\t' -f1 | tr '\n' ',' | sed 's/,$//g');
                echo -e "$fasta \t $NCBI \t $tax \t $motu \t $detect \t $map \t $percent \t $cog" >> classification.stats;
            done;
        done
        cd {config[path][root]}/{config[folder][abundance]}
        for folder in */;do 
            for file in $folder*.map;do 
                SAMP=$(paste $file | sed -n '1p' | cut -d ' ' -f1);
                BIN=$(paste $file | sed -n '3p' | cut -d ' ' -f1);
                LEN=$(paste $file | sed -n '12p' | cut -d ' ' -f4);
                MAG=$(paste $folder*_map.stats | sed -n '3p' | cut -d ' ' -f1);
                name=$(echo $file | sed 's/.map//g' | sed 's|^.*/||g');
                echo -n "$name ";
                awk -v samp="$SAMP" -v bin="$BIN" -v len="$LEN" 'BEGIN{printf bin/samp/len}';
                echo -n " "; 
                awk -v mag="$MAG" -v bin="$BIN" -v len="$LEN" 'BEGIN{print bin/mag/ len}';
            done;
        done > abundance.stats

        """


rule extractProteinBins:
    message:
        "Extract ORF annotated protein fasta files for each bin from reassembly checkm files."
    shell:
        """
        cd {config[path][root]}
        for folder in reassembled_bins/*/;do 
        for bin in $folder*reassembled_bins.checkm/bins/*;do 
        var=$(echo $bin/genes.faa | sed 's|reassembled_bins/||g'|sed 's|/reassembled_bins.checkm/bins||'|sed 's|/genes||g'|sed 's|/|_|g'|sed 's/permissive/p/g'|sed 's/orig/o/g'|sed 's/strict/s/g');
        cp $bin/*.faa /home/zorrilla/workspace/european/final_bins/$var;
        done;
        done
        """

rule carveme:
    input:
        bin = f'{config["path"]["root"]}/final_bins/{{binIDs}}.faa',
        media = config["dbs"]["carveme"]
    output:
        f'{config["path"]["root"]}/{config["folder"]["GEMs"]}/{{binIDs}}.xml'
    benchmark:
        f'{config["path"]["root"]}/benchmarks/{{binIDs}}.carveme.benchmark.txt'
    message:
        """
        Make sure that the input files are ORF annotated and preferably protein fasta.
        Will not work properly with raw genome fasta files.
        """
    shell:
        """
        set +u;source activate {config[envs][metabagpipes]};set -u
        mkdir -p $(dirname {output})
        cp {input.bin} {input.media} $TMPDIR
        cd $TMPDIR
        
        carve -g {config[params][carveMedia]} -v \
            --mediadb $(basename {input.media}) --fbc2 \
            -o $(echo $(basename {input.bin}) | sed 's/.faa/.xml/g') $(basename {input.bin})
        
        [ -f *.xml ] && mv *.xml $(dirname {output})
        """


rule modelVis:
    input:
        f'{config["path"]["root"]}/{config["folder"]["GEMs"]}'
    shell:
        """
        set +u;source activate memotenv;set -u;
        cd {input}
        for folder in ERR*;do 
        for model in $folder/*.xml;do 
        id=$(echo $model|sed 's|^.*/||g'|sed 's/.xml//g'); 
        mets=$(less $model| grep "species id="|cut -d ' ' -f 8|sed 's/..$//g'|sort|uniq|wc -l);
        rxns=$(less $model|grep -c 'reaction id=');
        genes=$(less $model|grep -c 'fbc:geneProduct fbc:id=');
        echo "$id $mets $rxns $genes" >> GEMs.stats;
        done;
        done
        
        Rscript {config[scripts][modelVis]}
        """


rule organizeGEMs:
    message:
        """
        Organizes GEMs into sample specific subfolders. 
        Necessary to run smetana per sample using the {IDs} wildcard.
        One liner (run from refined bins folder): 
        for folder in */;do mkdir -p ../GEMs/$folder;mv ../GEMs/$(echo $folder|sed 's|/||')_*.xml ../GEMs/$folder;done
        """
    shell:
        """
        cd {config[path][refined]}
        for folder in */;do 
        mkdir -p ../{config[path][GEMs]}; 
        mv ../{config[path][GEMs]}/$(echo $folder|sed 's|/||')_*.xml ../{config[path][GEMs]}/$folder;
        done
        """


rule smetana:
    input:
        f'{config["path"]["root"]}/{config["folder"]["GEMs"]}/{{IDs}}'
    output:
        f'{config["path"]["root"]}/{config["folder"]["SMETANA"]}/{{IDs}}_detailed.tsv'
    benchmark:
        f'{config["path"]["root"]}/benchmarks/{{IDs}}.smetana.benchmark.txt'
    shell:
        """
        set +u;source activate {config[envs][metabagpipes]};set -u
        mkdir -p {config[path][root]}/{config[folder][SMETANA]}
        cp {config[dbs][carveme]} {input}/*.xml $TMPDIR
        cd $TMPDIR
        
        smetana -o $(basename {input}) --flavor fbc2 \
            --mediadb media_db.tsv -m {config[params][smetanaMedia]} \
            --detailed \
            --solver {config[params][smetanaSolver]} -v *.xml
        
        cp *.tsv {config[path][root]} #safety measure for backup of results in case rule fails for some reason
        mv *.tsv $(dirname {output})
        """


rule interactionVis:
    input:
        f'{config["path"]["root"]}/{config["folder"]["SMETANA"]}'
    shell:
        """
        cd {input}
        mv media_db.tsv ../scripts/
        cat *.tsv|sed '/community/d' > smetana.all
        less smetana.all |cut -f2|sort|uniq > media.txt
        ll|grep tsv|awk '{print $NF}'|sed 's/_.*$//g'>samples.txt
        while read sample;do echo -n "$sample ";while read media;do var=$(less smetana.all|grep $sample|grep -c $media); echo -n "$var " ;done < media.txt; echo "";done < samples.txt > sampleMedia.stats
        """


rule memote:
    input:
        f'{config["path"]["root"]}/{config["folder"]["GEMs"]}/{{IDs}}'
    output:
        directory(f'{config["path"]["root"]}/{config["folder"]["memote"]}/{{IDs}}')
    benchmark:
        f'{config["path"]["root"]}/benchmarks/{{IDs}}.memote.benchmark.txt'
    shell:
        """
        set +u;source activate {config[envs][metabagpipes]};set -u
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
        bins = f'{config["path"]["root"]}/{config["folder"]["reassembled"]}/{{IDs}}/reassembled_bins',
        R1 = rules.qfilter.output.R1, 
        R2 = rules.qfilter.output.R2
    output:
        directory(f'{config["path"]["root"]}/{config["folder"]["GRiD"]}/{{IDs}}')
    benchmark:
        f'{config["path"]["root"]}/benchmarks/{{IDs}}.grid.benchmark.txt'
    shell:
        """
        set +u;source activate {config[envs][metabagpipes]};set -u
        cp -r {input.bins} {input.R1} {input.R2} $TMPDIR
        cd $TMPDIR
        cat *.gz > $(basename $(dirname {input.bins})).fastq.gz
        rm $(basename {input.R1}) $(basename {input.R2})
        mkdir MAGdb out
        update_database -d MAGdb -g $(basename {input.bins}) -p MAGdb
     	rm -r $(basename {input.bins})
        grid multiplex -r . -e fastq.gz -d MAGdb -p -c 0.2 -o out -n {config[cores][grid]}
        rm $(basename $(dirname {input.bins})).fastq.gz
        mkdir {output}
        mv out/* {output}
        """
