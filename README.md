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

## Installation

```bash
git clone https://github.com/franciscozorrilla/metaGEM.git
```

### Conda

#### metaGEM

A [conda](https://conda.io/en/latest/) specification file is provided to install required packages inside an
environment called `metaGEM`.

```bash
conda env create -f metaGEM_env.yml
source activate metaGEM
```

#### metaWRAP

Bin refinement and bin reassembly make use of metaWRAP modules. To avoid package conflicts, set up metaWRAP in its own environment using the provided conda specification file.

```bash
conda env create -f metaWRAP_env.yml
source activate metawrap
```

CheckM is used extensively to evaluate the output of various itntermediate steps. Although the CheckM package is installed in the `metawrap` environment, the user is required to download the CheckM database and run `checkm data setRoot <db_dir>` as outlined in the [CheckM installation guide](https://github.com/Ecogenomics/CheckM/wiki/Installation#how-to-install-checkm).

##### CPLEX

GEM reconstruction and GEM community simulations require the IBM CPLEX solver, which is [free to download with an academic license](https://developer.ibm.com/docloud/blog/2019/07/04/cplex-optimization-studio-for-students-and-academics/). Refer to the [CarveMe](https://carveme.readthedocs.io/en/latest/installation.html) and [SMETANA](https://smetana.readthedocs.io/en/latest/installation.html) installation instructions for further information or troubleshooting. Note: CPLEX v.12.8 is recommended.

### Singularity

Alternatively, metaGEM can be installed with the provided [Singularity](https://sylabs.io/docs/) recipe file.

```bash
sudo singularity --verbose build metaGEM.simg Singularity
```
