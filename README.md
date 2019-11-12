# metaBAGpipes
### A Snakemake-based workflow to generate high quality MAGs, reconstruct GEMs, and perform community metabolic interaction simulations on HPC cluster

## Use cases:

1. Assembly using [metaSPAdes](https://github.com/ablab/spades).
2. Binning using:
      * [CONCOCT](https://github.com/BinPro/CONCOCT).
      * [metabat2](https://bitbucket.org/berkeleylab/metabat/src/master/).
      * [maxbin2](https://sourceforge.net/projects/maxbin2/).
3. Bin refinement using [metaWRAP](https://github.com/bxlab/metaWRAP):
      * reconciles binning output from different tools using [CheckM](https://github.com/Ecogenomics/CheckM) and [binning_refiner](https://github.com/songweizhi/Binning_refiner).
4. Bin reassembly using [metaWRAP](https://github.com/bxlab/metaWRAP):
      * extracts reads corresponding to bins from each sample and re-assembles using [SPAdes](https://github.com/ablab/spades).
5. MAG classification using [classify genomes](https://github.com/AlessioMilanese/classify-genomes) based on mOTUs2.
6. MAG abundance:
   * using [mOTUs2](https://github.com/motu-tool/mOTUs_v2).
   * based on MAG-sample mapping using [bwa](https://github.com/lh3/bwa).
7. GEM reconstruction using [CarveMe](https://github.com/cdanielmachado/carveme).
8. GEM QC using [memote](https://github.com/opencobra/memote).
9. GEM community simulations using [SMETANA](https://github.com/cdanielmachado/smetana).
10. MAG growth rate estimates using [GRiD](https://github.com/ohlab/GRiD).

## Toy dataset

As a tutorial, and to verify that your metaBAGpipes installation is working correctly, we provide a [toy dataset](https://zenodo.org/record/3534949#.XclQriV7lTZ). Expected output can be found in the Tutorial folder.

## Installation

```
git clone https://github.com/franciscozorrilla/metaBAGpipes.git
```

### Conda environments

#### metaBAGpipes

A [conda](https://conda.io/en/latest/) specification file is provided to install required packages inside an
environment called `metabagpipes`.

```
conda env create -f metaBAGpipes_env.yml
source activate metabagpipes
```

##### CPLEX

GEM reconstruction (step 7) and GEM community simulations (step 9) require the IBM CPLEX solver, which is [free to download with an academic license](https://developer.ibm.com/docloud/blog/2019/07/04/cplex-optimization-studio-for-students-and-academics/). Refer to the [CarveMe](https://carveme.readthedocs.io/en/latest/installation.html) and [SMETANA](https://smetana.readthedocs.io/en/latest/installation.html) installation instructions for further information or troubleshooting.

#### metaWRAP

Bin refinement (step 3) and bin reassembly (step 4) make use of metaWRAP modules. To avoid package conflicts, set up metaWRAP in its own environment using the provided conda specification file.

```
conda env create -f metaWRAP_env.yml
source activate metawrap
```

##### CheckM

CheckM is used extensively to evaluate the output of various itntermediate steps. Although the CheckM package is installed in the `metawrap` environment, the user is required to download the CheckM database and run `checkm data setRoot <db_dir>` as outlined in the [CheckM installation guide](https://github.com/Ecogenomics/CheckM/wiki/Installation#how-to-install-checkm).

### Singularity

A [Singularity](https://sylabs.io/docs/) recipe files is provided to build an image that can be used with HPC clusters.

```
sudo singularity --verbose build metabagpipes.simg Singularity
```

## Tutorial

### 0. Setup

Create tutorial directory and move `Snakefile`, `config.yaml`, and `cluster_config.json` to working directory:

```
mkdir -p tutorial
mv Snakefile config.yaml cluster_config.json tutorial
cd tutorial
pwd
```

#### Configure config.yaml

The `config.yaml` contains the absolute path of your working directory, folder names, absolute paths to scripts/databases, and cores/parameters used by individual rules in the Snakefile. Most importantly, replace the root path (line 2) with the absolute path of your working directory in the `config.yaml` file.

```
path:
    root: /home/zorrilla/workspace/githubReadmeTutorial
folder:
    data: dataset
    logs: logs
    assemblies: assemblies
    scripts: scripts
    concoctInput: concoct_input
    concoctOutput: concoct_output
    maxbin: maxbin_output
    metabat: metabat_output
    refined: refined_bins
    reassembled: reassembled_bins
    classification: classification
    abundance: abundance
    GRiD: GRiD
    GEMs: GEMs
    SMETANA: SMETANA
    memote: memote
scripts:
    kallisto2concoct: /home/zorrilla/workspace/tutorial/scripts/kallisto2concoct.py
    binFilter: /home/zorrilla/workspace/tutorial/scripts/binFilter.py
dbs:
    carveme: /home/zorrilla/workspace/tutorial/scripts/media_db.tsv
    toy: /home/zorrilla/workspace/tutorial/scripts/download_toydata.txt
cores:
    metaspades: 24
    kallisto: 24
    concoct: 24
    metabat: 24
    maxbin: 24
    refine: 24
    reassemble: 24
    classify: 2
    abundance: 16
    carveme: 4
    smetana: 12
    memote: 4
    grid: 24
params:
    cutfasta: 10000
    concoct: 800
    refineMem: 1600
    refineComp: 50
    refineCont: 10
    reassembleMem: 1600
    reassembleComp: 50
    reassembleCont: 10
    carveMedia: M8
    smetanaMedia: M1,M2,M3,M4,M5,M7,M8,M9,M10,M11,M13,M14,M15A,M15B,M16
    smetanaSolver: CPLEX
envs:
    metabagpipes: metabagpipes
    metawrap: metawrap
```

The `folders` section should be left as is. However, ensure that the `scripts` and `dbs` point to the correct folders/files. The `cores` section can be modified to best suit the architecture of your cluster. Note that these parameters are only used by the `Snakefile`, so the `cluster_config.json` file also needs to be modified for each batch of jobs.

#### Configure Snakefile wildcards

To test that your snakemake installation and Snakefile are working properly, run the `createFolders` snakemake rule. These folders can alternatively be generated later on during the execution of each individual rule.

```
snakemake createFolders
```

Next, download the [toy dataset](https://zenodo.org/record/3534949#.XclQriV7lTZ) into the `dataset` folder. This can be done manually, or more conveniently, using the `downloadToy` snakemake rule.

```
snakemake downloadToy
```

Organize reads into sample specific sub-directories. This is required as metaBAGpipes uses snakemake wildcards based on these subfolders.

```
snakemake organizeData
```

Ensure that the line 6 of your Snakefile is indeed pointing to the dataset folder with sample specific subfolders. In the case of this tutorial:

```
IDs = sorted([os.path.splitext(val)[0] for val in (glob.glob('dataset/*'))])
```

#### Configure cluster_config.json

The main body of the cluster_config.json file should look like this:

```
{
"__default__" : {
        "account" : "patil",
        "time" : "7-00:00:00",
        "n" : 16,
        "tasks" : 1,
        "mem" : 60G,
        "cpusPerTask" : 16,
        "name"      : "DL.{rule}",
        "output"    : "logs/{wildcards}.%N.{rule}.out.log",
},
}
```

Configure the `cluster_config.json` file by editing the account field (line 3). The majority of the pipeline will be run through the cluster using the last line in the `cluster_config.json` file:

```
nohup snakemake all -j 200 -k --cluster-config cluster_config.json -c "sbatch -A {cluster.account} -t {cluster.time} -n {cluster.n} --ntasks {cluster.tasks} --cpus-per-task {cluster.cpusPerTask} --output {cluster.output}" &
```

Note that the `n` (line 5) and `cpusPerTask` (line 8) should always match. Read the [slurm documentation](https://slurm.schedmd.com/documentation.html) to learn about different possible flags.

#### Note: rule all

Since snakemake rules with wildcards cannot be target rules, we use the input of `rule all` to expand our wildcards. For example, to run SMETANA one would do:

```
rule all:
    input:
        expand(config["path"]["root"]+"/"+config["folder"]["SMETANA"]+"/{IDs}_detailed.tsv", IDs = IDs)
    shell:
        """
        echo {input}
        """
```

### 1. Assembly

We use the metaSPAdes assembler which contains its own internal quality control module. If you wish to use another assembler such as [megahit](https://github.com/voutcn/megahit), we recommend pre-processing using [fastp](https://github.com/OpenGene/fastp) or [trimmomatic](https://github.com/timflutre/trimmomatic).

To run the assembly step, copy the output of `rule metaspades`, and insert it into the input for `rule all`:

```
rule all:
    input:
        expand(config["path"]["root"]+"/"+config["folder"]["assemblies"]+"/{IDs}/contigs.fasta.gz", IDs = IDs)
    shell:
        """
        echo {input}
        """
```

Next, edit the `n` and `cpusPerTask` fields of the `cluster_config.json` file to specify the desired number of cores for the assembly jobs.

Finally, run the following snippet of code to submit our assembly jobs to the cluster scheduler. Note that the `-j` flag specifies the number of parallel jobs to be submitted. In the toy dataset we have 3 samples, therefore we specify:

```
nohup snakemake all -j 3 -k --cluster-config cluster_config.json -c "sbatch -A {cluster.account} -t {cluster.time} -n {cluster.n} --ntasks {cluster.tasks} --cpus-per-task {cluster.cpusPerTask} --output {cluster.output}" &
```

The snakemake output should be stored to `nohup.out`, while the output for the individual jobs can be found in the `logs` folder. Furthermore, benchmarks containing job-specific metrics are generated in the `benchmarks` folder.

One can easily check the status of running/pending jobs by running:

```
squeue -u <USER>
```

### 2. Binning

Run CONCOT, metabat2, and maxbin2 individually. The output of these different tools will later be reconciled using the metaWRAP bin refinement module.

#### CONCOCT

CONCOCT requires coverage information across samples. To generate this input, we use `rule kallisto` to map each sample to each assembly.

Edit the input of `rule all` to expand the output of `rule kallisto`:

```
rule all:
    input:
        expand(config["path"]["root"]+"/"+config["folder"]["concoctInput"]+"/{IDs}_concoct_inputtableR.tsv", IDs = IDs)
    shell:
        """
        echo {input}
        """
```

Edit `cluster_config.json` if necessary and run snakemake using:

```
nohup snakemake all -j 3 -k --cluster-config cluster_config.json -c "sbatch -A {cluster.account} -t {cluster.time} -n {cluster.n} --ntasks {cluster.tasks} --cpus-per-task {cluster.cpusPerTask} --output {cluster.output}" &
```

Once these jobs finish, edit the input of `rule all` to expand the output of `rule concoct`:

```
rule all:
    input:
        expand(config["path"]["root"]+"/"+config["folder"]["concoctOutput"]+"/{IDs}/{IDs}.concoct-bins", IDs = IDs)
    shell:
        """
        echo {input}
        """
```

Note that we omit the `directory()` flag from the `rule concoct` output when we expand wildcards in `rule all`, as it is not necessary for rule inputs.

```
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
        set +u;source activate {config[envs][metabagpipes]};set -u;
        mkdir -p $(dirname {output})
        mkdir -p $(dirname $(dirname {output}))
        cp {input.contigs} {input.table} $TMPDIR
        cd $TMPDIR
        gunzip $(basename {input.contigs})
        cut_up_fasta.py -c {config[params][cutfasta]} -o 0 -m contigs.fasta > metaspades_c10k.fa
        concoct --coverage_file {input.table} --composition_file metaspades_c10k.fa -b $(basename $(dirname {output})) -t {config[cores][concoct]} -c {config[params][concoct]}
        merge_cutup_clustering.py $(basename $(dirname {output}))_clustering_gt1000.csv > $(basename $(dirname {output}))_clustering_merged.csv
        mkdir -p $(basename {output})
        extract_fasta_bins.py contigs.fasta $(basename $(dirname {output}))_clustering_merged.csv --output_path $(basename {output})
        checkm lineage_wf -x fa -t {config[cores][concoct]} --pplacer_threads {config[cores][concoct]} --file $(basename $(dirname {output})).tab --tab_table $(basename {output}) .
        mv $(basename {output}) *.log *.txt *.csv *.tab $(dirname {output})
        """
```        

Make sure to similarly remove the `directory()` flag when expanding wildcards in `rule all` for all future rules.

Edit `cluster_config.json` if necessary and run snakemake using:

```
nohup snakemake all -j 3 -k --cluster-config cluster_config.json -c "sbatch -A {cluster.account} -t {cluster.time} -n {cluster.n} --ntasks {cluster.tasks} --cpus-per-task {cluster.cpusPerTask} --output {cluster.output}" &
```

#### metabat2

Edit the input of `rule all` to expand the output of `rule metabat`, edit `cluster_config.json` if necessary, and run snakemake just as before.

#### maxbin2

Edit the input of `rule all` to expand the output of `rule maxbin`, edit `cluster_config.json` if necessary, and run snakemake.

### 3. Bin refinement

Edit the input of `rule all` to expand the output of `rule binRefine`, edit `cluster_config.json` if necessary, and run snakemake. Note that this step uses the `metawrap` conda environment instead of the `metabagpipes` environment. Ensure that CheckM is correctly installed and the downloaded database is configured as described in the Installation section.

### 4. Bin reassembly

Edit the input of `rule all` to expand the output of `rule binReassemble`, edit `cluster_config.json` if necessary, and run snakemake. Note that this step uses the `metawrap` conda environment instead of the `metabagpipes` environment. Again, ensure that CheckM is correctly installed and the downloaded database is configured as described in the Installation section.

### 5. MAG classification

Edit the input of `rule all` to expand the output of `rule classifyGenomes`, edit `cluster_config.json` if necessary, and run snakemake.

### 6. MAG abundance

Edit the input of `rule all` to expand the output of `rule abundance`, edit `cluster_config.json` if necessary, and run snakemake.

### 7. GEM reconstruction

### 8. GEM quality control

### 9. Community simulations

### 10. GEM/MAG growth rate

## Abstract
metaBAGpipes integrates an array of existing bioinformatics and metabolic modeling tools using Snakemake, for the purpose of interrogating social interactions in bacterial communities of the human gut microbiome. From WGS metagenomic datasets, metagenome assembled genomes (MAGs) are reconstructed, which are then converted into genome-scale metabolic models (GEMs) for *in silico* simulations of cross feeding interactions within sample based communities. Abundance estimates for community members are estimated by mapping metagenomic samples to the generated MAGs, which are used in combination with the simulated cross feeding interactions for the generation of explanatory and statistically significant linear models. We conclude that there is indeed a correlation, ranging from weak to moderate, between gut microbiome members’ abundance and set of metabolic cross-feeding interactions across samples. A more comprehensive analysis incorporating multiple datasets needs to be conducted to strengthen and expand the findings of this work.

## Significance

To our best knowledge, the work presented in this thesis project represents the first application of sample-specific gut microbiome community FBA simulations using MAG-based GEMs. While there are limitations to this approach, as discussed in the thesis write-up, it is conceivable that a fine-tuned version of such a pipeline could be used to study, develop, and evualuate personalized medicine treatments. 137 WGS gut metagenome samples were processed in this study, a relatively modest amount of data. By analyzing a greater number of samples, further machine learning methods can be employed to interrogate the network architechture of microbial communities and facilitiate data driven hypothesis generation.

Dataset used:
  * Karlsson, Fredrik H., et al. “Gut Metagenome in European Women with Normal, Impaired and Diabetic Glucose Control.” *Nature*, vol.498, no.7452, 2013, pp.99–103. , doi:10.1038/nature12198.

This repository is administered by Francisco Zorrilla ([@franciscozorrilla](https://github.com/franciscozorrilla/)), Structural and Computational Biology Unit, EMBL. metaBAGpipes was developed throughout my Master's thesis project at the Systems and Synthetic Biology division of Chalmers Univeristy of Technology, under the supervision of Aleksej Zelezniak.

  * Last update: 11-11-2019
