# metaGEM

> A Snakemake-based workflow to generate high quality MAGs, reconstruct GEMs, and perform community metabolic interaction simulations on HPC clusters.

![metawrapfigs_v2 002](https://user-images.githubusercontent.com/35606471/103545679-ceb71580-4e99-11eb-9862-084115121980.jpeg)

metaGEM integrates an array of existing bioinformatics and metabolic modeling tools using Snakemake, for the purpose of interrogating social interactions in bacterial communities of the human gut microbiome. From WMGS datasets, metagenome assembled genomes (MAGs) are reconstructed, which are then converted into genome-scale metabolic models (GEMs) for *in silico* simulations of cross feeding interactions within sample based communities. Additional outputs include abundance estimates and taxonomic annotations.

## Use cases:

1. Raw read QC/trimming using [fastp](https://github.com/OpenGene/fastp).
2. Assembly using [megahit](https://github.com/voutcn/megahit).
3. Binning using:
      * [CONCOCT](https://github.com/BinPro/CONCOCT).
      * [metabat2](https://bitbucket.org/berkeleylab/metabat/src/master/).
      * [maxbin2](https://sourceforge.net/projects/maxbin2/).
4. Bin refinement using [metaWRAP](https://github.com/bxlab/metaWRAP):
      * reconciles binning output from different tools using [CheckM](https://github.com/Ecogenomics/CheckM) and [binning_refiner](https://github.com/songweizhi/Binning_refiner).
5. Bin reassembly using [metaWRAP](https://github.com/bxlab/metaWRAP):
      * extracts reads corresponding to bins from each sample and re-assembles using [SPAdes](https://github.com/ablab/spades).
6. MAG classification using [classify genomes](https://github.com/AlessioMilanese/classify-genomes) based on mOTUs2.
7. MAG abundance:
   * using [mOTUs2](https://github.com/motu-tool/mOTUs_v2).
   * based on MAG-sample mapping using [bwa](https://github.com/lh3/bwa).
8. GEM reconstruction using [CarveMe](https://github.com/cdanielmachado/carveme).
9. GEM QC using [memote](https://github.com/opencobra/memote).
10. GEM community simulations using [SMETANA](https://github.com/cdanielmachado/smetana).

Experimental:

11. MAG growth rate estimates using [GRiD](https://github.com/ohlab/GRiD).

### Usage

bash metaGEM.sh [-t TASK] [-j NUMBER OF JOBS] [-c NUMBER OF CORES]

```

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
                            kallisto
                            concoct
                            metabat
                            maxbin
                            binRefine
                            binReassemble
                            classifyGenomes
                            abundance
                            moveBins
                            carveme
                            organizeGems
                            smetana
                            memote   
                            grid
                        VISUALIZATION (in development)
                            assemblyVis
                            binningVis
                            taxonomyVis
                            abundanceVis
                            modelVis
                            interactionVis
                            growthVis
  -j, --nJobs       Specify number of jobs to run in parallel
  -c, --nCores      Specify number of cores per job
  -h, --help        Display this help and exit

Example: bash metaGEM.sh -t createFolders -j 1 -c 1

```

## Installation

```bash
git clone https://github.com/franciscozorrilla/metaGEM.git
```

### Conda environments

#### metaGEM

A [conda](https://conda.io/en/latest/) specification file is provided to install required packages inside an
environment called `metaGEM`.

```bash
conda env create -f metaGEM_env.yml
source activate metaGEM
```

##### CPLEX

GEM reconstruction (step 7) and GEM community simulations (step 9) require the IBM CPLEX solver, which is [free to download with an academic license](https://developer.ibm.com/docloud/blog/2019/07/04/cplex-optimization-studio-for-students-and-academics/). Refer to the [CarveMe](https://carveme.readthedocs.io/en/latest/installation.html) and [SMETANA](https://smetana.readthedocs.io/en/latest/installation.html) installation instructions for further information or troubleshooting. Note: CPLEX v.12.8 is recommended.

#### metaWRAP

Bin refinement (step 3) and bin reassembly (step 4) make use of metaWRAP modules. To avoid package conflicts, set up metaWRAP in its own environment using the provided conda specification file.

```bash
conda env create -f metaWRAP_env.yml
source activate metawrap
```

##### CheckM

CheckM is used extensively to evaluate the output of various itntermediate steps. Although the CheckM package is installed in the `metawrap` environment, the user is required to download the CheckM database and run `checkm data setRoot <db_dir>` as outlined in the [CheckM installation guide](https://github.com/Ecogenomics/CheckM/wiki/Installation#how-to-install-checkm).

### Singularity

A [Singularity](https://sylabs.io/docs/) recipe files is provided to build an image that can be used with HPC clusters.

```bash
sudo singularity --verbose build metaGEM.simg Singularity
```
