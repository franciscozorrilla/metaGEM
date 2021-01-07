# metaGEM

> A Snakemake-based workflow to generate high quality MAGs, reconstruct GEMs, and perform community metabolic interaction simulations on HPC clusters.

![metawrapfigs_v2 002](https://user-images.githubusercontent.com/35606471/103545679-ceb71580-4e99-11eb-9862-084115121980.jpeg)

metaGEM integrates an array of existing bioinformatics and metabolic modeling tools using Snakemake, for the purpose of interrogating social interactions in bacterial communities of the human gut microbiome. From WMGS datasets, metagenome assembled genomes (MAGs) are reconstructed, which are then converted into genome-scale metabolic models (GEMs) for *in silico* simulations of cross feeding interactions within sample based communities. Additional outputs include abundance estimates and taxonomic annotations.

## Workflow

0. metaGEM setup
1. Quality filter reads with [fastp](https://github.com/OpenGene/fastp)
2. Assembly with [megahit](https://github.com/voutcn/megahit)
3. Draft bin sets with [CONCOCT](https://github.com/BinPro/CONCOCT),[MaxBin2](https://sourceforge.net/projects/maxbin2/), and [MetaBAT2](https://sourceforge.net/projects/maxbin2/)
4. Refine & reassemble bins with [metaWRAP](https://github.com/bxlab/metaWRAP)
5. Taxonomic assignment with [GTDB-tk](https://github.com/Ecogenomics/GTDBTk)
6. Relative abundances with [bwa](https://github.com/lh3/bwa) and [samtools](https://github.com/samtools/samtools)
7. Reconstruct & evaluate genome-scale metabolic models with [CarveMe](https://github.com/cdanielmachado/carveme) and [memote](https://github.com/opencobra/memote)
8. Species metabolic coupling analysis with [SMETANA](https://github.com/cdanielmachado/smetana)
9. Growth rate estimation with [GRiD](https://github.com/ohlab/GRiD)
10. Pangenome analysis with [roary](https://github.com/sanger-pathogens/Roary)
11. Eukaryotic draft bins with [EukRep](https://github.com/patrickwest/EukRep) and [EukCC](https://github.com/Finn-Lab/EukCC)

### Usage

```
_________________________________________________________________________/\\\\\\\\\\\\___/\\\\\\\\\\\\\\\___/\\\\____________/\\\\_        
 _______________________________________________________________________/\\\//////////___\/\\\///////////___\/\\\\\\________/\\\\\\_       
  __________________________________________/\\\________________________/\\\______________\/\\\______________\/\\\//\\\____/\\\//\\\_      
   ____/\\\\\__/\\\\\________/\\\\\\\\____/\\\\\\\\\\\___/\\\\\\\\\_____\/\\\____/\\\\\\\__\/\\\\\\\\\\\______\/\\\\///\\\/\\\/_\/\\\_     
    __/\\\///\\\\\///\\\____/\\\/////\\\__\////\\\////___\////////\\\____\/\\\___\/////\\\__\/\\\///////_______\/\\\__\///\\\/___\/\\\_    
     _\/\\\_\//\\\__\/\\\___/\\\\\\\\\\\______\/\\\_________/\\\\\\\\\\___\/\\\_______\/\\\__\/\\\______________\/\\\____\///_____\/\\\_   
      _\/\\\__\/\\\__\/\\\__\//\\///////_______\/\\\_/\\____/\\\/////\\\___\/\\\_______\/\\\__\/\\\______________\/\\\_____________\/\\\_  
       _\/\\\__\/\\\__\/\\\___\//\\\\\\\\\\_____\//\\\\\____\//\\\\\\\\/\\__\//\\\\\\\\\\\\/___\/\\\\\\\\\\\\\\\__\/\\\_____________\/\\\_ 
        _\///___\///___\///_____\//////////_______\/////______\////////\//____\////////////_____\///////////////___\///______________\///__
        
        
Usage: bash metaGEM.sh [-t|--task TASK] 
                       [-j|--nJobs NUMBER OF JOBS] 
                       [-c|--cores NUMBER OF CORES] 
                       [-m|--mem GB RAM] 
                       [-h|--hours MAX RUNTIME]

Snakefile wrapper/parser for metaGEM. 

Options:
  -t, --task        Specify task to complete:

                        SETUP
                            createFolders
                            downloadToy
                            organizeData

                        WORKFLOW
                            fastp 
                            megahit 
                            crossMap 
                            concoct 
                            metabat
                            maxbin 
                            binRefine 
                            binReassemble 
                            extractProteinBins
                            carveme
                            memote
                            organizeGEMs
                            smetana
                            extractDnaBins
                            gtdbtk
                            abundance 
                            grid
                            prokka
                            roary

                        VISUALIZATION (in development)
                            qfilterVis
                            assemblyVis
                            binningVis
                            taxonomyVis
                            modelVis
                            interactionVis
                            growthVis

  -j, --nJobs       Specify number of jobs to run in parallel
  -c, --nCores      Specify number of cores per job
  -m, --mem         Specify memory in GB required for job
  -h, --hours       Specify number of hours to allocated to job runtime
```

## Automated installation

Clone this repository to your HPC or local computer and run the `env_setup.sh` script:

```
git clone https://github.com/franciscozorrilla/metaGEM.git
cd metaGEM
bash env_setup.sh
```

This script will set up 3 conda environments, `metagem`, `metawrap`, and `prokkaroary`, which will be activated as required by Snakemake jobs.

### CheckM

CheckM is used extensively to evaluate the output of various intermediate steps. Although the CheckM package is installed in the `metawrap` environment, the user is required to download the CheckM database and run `checkm data setRoot <db_dir>` as outlined in the [CheckM installation guide](https://github.com/Ecogenomics/CheckM/wiki/Installation#how-to-install-checkm).

### CPLEX

Unfortunately CPLEX cannot be automatically installed in the `env_setup.sh` script, you must install this dependency manually within the metagem conda environment. GEM reconstruction and GEM community simulations require the IBM CPLEX solver, which is [free to download with an academic license](https://developer.ibm.com/docloud/blog/2019/07/04/cplex-optimization-studio-for-students-and-academics/). Refer to the [CarveMe](https://carveme.readthedocs.io/en/latest/installation.html) and [SMETANA](https://smetana.readthedocs.io/en/latest/installation.html) installation instructions for further information or troubleshooting. Note: CPLEX v.12.8 is recommended.

## Manual installation

You can manually set up the environments with the following chunks of code.

### metaGEM

```
conda create -n metagem mamba
source activate metagem
mamba install python snakemake fastp megahit bwa samtools=1.9 kallisto concoct=1.1 metabat2 maxbin2 gtdbtk eukrep eukcc smeg motus
pip install --user memote carveme smetana
```

### metaWRAP

```
conda create -n metawrap
source activate metawrap
conda install -c ursky metawrap-mg=1.3.2
```

### prokka-roary

```
conda create -n prokkaroary
source activate prokkaroary
conda install prokka roary
```

## Singularity

Alternatively, metaGEM can be installed with the provided [Singularity](https://sylabs.io/docs/) recipe file.

```bash
sudo singularity --verbose build metaGEM.simg Singularity
```

## Resources

### Unseen bio demo

metaGEM can be used to explore your own gut microbiome based on at-home-test-kit seqencing data from services such as [unseen bio](https://unseenbio.com/). The [following demo](https://github.com/franciscozorrilla/unseenbio_metaGEM) showcases the metaGEM workflow on two unseenbio samples.

### Toy dataset tutorial

Although compiled while metaGEM was still in development, [this slightly outdated tutorial](https://github.com/franciscozorrilla/metaGEM/blob/master/Tutorial/tutorial.md) can provide additional insight as to how metaGEM is implemented and used. The [toy dataset](https://zenodo.org/record/3534949#.XclQriV7lTZ) was generated by subsetting real shotgun [sequencing samples](https://www.ncbi.nlm.nih.gov/sra?term=ERP002469) from the human gut using the [seqtk tool](https://github.com/lh3/seqtk):

```
seqtk sample -s100 sample_X.fastq.gz 3000000 > subset_X.fastq.gz
```

The idea was to design a small and portable dataset amenable to unit testing, however the size reduction came at the expense of performance. In retrospect, it may have been a better idea to [generate a synthetic/artificial metagenomic dataset](https://github.com/CAMI-challenge/CAMISIM) from a handful of genomes.

### Publications

The metaGEM workflow has been used in some capacity in the following publications:

```
Plastic-degrading potential across the global microbiome correlates with recent pollution trends
Jan Zrimec, Mariia Kokina, Sara Jonasson, Francisco Zorrilla, Aleksej Zelezniak
bioRxiv 2020.12.13.422558; doi: https://doi.org/10.1101/2020.12.13.422558 
```

## Please cite

```
metaGEM: reconstruction of genome scale metabolic models directly from metagenomes
Francisco Zorrilla, Kiran R. Patil, Aleksej Zelezniak
bioRxiv 2020.12.31.424982; doi: https://doi.org/10.1101/2020.12.31.424982 
```
