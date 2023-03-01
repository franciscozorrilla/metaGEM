# 💎 `metaGEM`

> **Note** 
> An easy-to-use workflow for generating context specific genome-scale metabolic models and predicting metabolic interactions within microbial communities directly from metagenomic data.

[![Nucleic Acids Research](https://img.shields.io/badge/Nucleic%20Acids%20Research-10.1093%2Fnar%2Fgkab815-critical)](https://academic.oup.com/nar/advance-article/doi/10.1093/nar/gkab815/6382386)
[![bioRxiv](https://img.shields.io/badge/bioRxiv-10.1101%2F2020.12.31.424982%20-B31B1B)](https://www.biorxiv.org/content/10.1101/2020.12.31.424982v2.full)
[![Build Status](https://app.travis-ci.com/franciscozorrilla/metaGEM.svg?branch=master)](https://app.travis-ci.com/github/franciscozorrilla/metaGEM)
[![GitHub license](https://img.shields.io/github/license/franciscozorrilla/metaGEM)](https://github.com/franciscozorrilla/metaGEM/blob/master/LICENSE)
[![Gitter chat](https://badges.gitter.im/gitterHQ/gitter.png)](https://gitter.im/metaGEM/community)
[![DOI](https://img.shields.io/badge/Zenodo-10.5281%2F4707723-blue)](https://zenodo.org/badge/latestdoi/137376259)
[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/drive/1I1S8AoGuJ9Oc2292vqAGTDmZcbnolbuj#scrollTo=awiAaVwSF5Fz)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/metagem/badges/version.svg)](https://anaconda.org/bioconda/metagem)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/metagem/badges/latest_release_date.svg)](https://anaconda.org/bioconda/metagem)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/metagem/badges/downloads.svg)](https://anaconda.org/bioconda/metagem)

![metawrapfigs_final4 001](https://user-images.githubusercontent.com/35606471/116543667-0d0f8f00-a8e6-11eb-835c-bc1fe935f43e.png)

`metaGEM` is a Snakemake workflow that integrates an array of existing bioinformatics and metabolic modeling tools, for the purpose of predicting metabolic interactions within bacterial communities of microbiomes. From whole metagenome shotgun datasets, metagenome assembled genomes (MAGs) are reconstructed, which are then converted into genome-scale metabolic models (GEMs) for *in silico* simulations. Additional outputs include abundance estimates, taxonomic assignment, growth rate estimation, pangenome analysis, and eukaryotic MAG identification.

## ⚙️ Installation

You can start using `metaGEM` on your cluster with just one line of code with the [mamba package manager](https://github.com/mamba-org/mamba)

```
mamba create -n metagem -c bioconda metagem
```

This will create an environment called `metagem` and start installing dependencies. Please consult the `config/README.md` page for more detailed setup instructions.

[![installation](https://img.shields.io/badge/metaGEM-config-%2331a354)](https://github.com/franciscozorrilla/metaGEM/tree/master/config)

## 🔧 Usage

Clone this repo

```
git clone https://github.com/franciscozorrilla/metaGEM.git && cd metaGEM/workflow
```

Run `metaGEM` without any arguments to see usage instructions:

```
bash metaGEM.sh
```
```
Usage: bash metaGEM.sh [-t|--task TASK] 
                       [-j|--nJobs NUMBER OF JOBS] 
                       [-c|--cores NUMBER OF CORES] 
                       [-m|--mem GB RAM] 
                       [-h|--hours MAX RUNTIME]
                       [-l|--local]
                       
 Options:
  -t, --task        Specify task to complete:

                        SETUP
                            createFolders
                            downloadToy
                            organizeData
                            check

                        CORE WORKFLOW
                            fastp 
                            megahit 
                            crossMapSeries
                            kallistoIndex
                            crossMapParallel
                            kallisto2concoct 
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

                        BONUS
                            grid
                            prokka
                            roary
                            eukrep
                            eukcc

                        VISUALIZATION (in development)
                            stats
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
  -l, --local       Run jobs on local machine for non-cluster usage
```

## 🧉 Try it now

You can set up and use `metaGEM` on the cloud by following along the google colab notebook. 

[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/drive/1I1S8AoGuJ9Oc2292vqAGTDmZcbnolbuj#scrollTo=awiAaVwSF5Fz)

Please note that google colab does not provide the computational resources necessary to fully run `metaGEM` on a real dataset. This notebook demonstrates how to set up and use `metaGEM` by perfoming the first steps in the workflow on a toy dataset.

## 💩 Tutorials

`metaGEM` can be used to explore your own gut microbiome sequencing data from at-home-test-kit services such as [unseen bio](https://unseenbio.com/). The following tutorial showcases the `metaGEM` workflow on two unseenbio samples.

[![Tutorial](https://img.shields.io/badge/metaGEM-Tutorial-%23d8b365)](https://github.com/franciscozorrilla/unseenbio_metaGEM)

For an introductory metabolic modeling tutorial, refer to the resources compiled for the [EMBOMicroCom: Metabolite and species dynamics in microbial communities](https://www.embl.org/about/info/course-and-conference-office/events/mcd22-01/) workshop in 2022.

[![Tutorial3](https://img.shields.io/badge/MicroCom-Tutorial-green)](https://github.com/franciscozorrilla/EMBOMicroCom)

For a more advanced tutorial, check out the resources we put together for the [SymbNET: from metagenomics to metabolic interactions](https://www.ebi.ac.uk/training/events/symbnet-2022/) course in 2022.

[![Tutorial2](https://img.shields.io/badge/SymbNET-Tutorial-red)](https://github.com/franciscozorrilla/SymbNET)

## 🏛️ Wiki

Refer to the wiki for additional usage tips, frequently asked questions, and implementation details.

[![wiki](https://img.shields.io/badge/metaGEM-Wiki-blue)](https://github.com/franciscozorrilla/metaGEM/wiki)

## 📦 Datasets

* You can access the metaGEM-generated results for the publication [here](https://github.com/franciscozorrilla/metaGEM_paper).
```
    🧪 Small communities of gut microbes from lab cultures
    💩 Real gut microbiome samples from Swedish diabetes paper
    🪴 Plant-associated soil samples from Chinese rhizobiome study
    🌏 Bulk-soil samples from Australian biodiversity analysis
    🌊 Ocean water samples from global TARA Oceans expeditions
```
* Additionally, you can access metaGEM-generated results from a reanalysis of recently published ancient metagenomes [here](https://zenodo.org/record/7414438#.Y5HSFYLP3bs).

## 🐍 Workflow

### Core

1. Quality filter reads with [fastp](https://github.com/OpenGene/fastp)
2. Assembly with [megahit](https://github.com/voutcn/megahit)
3. Draft bin sets with [CONCOCT](https://github.com/BinPro/CONCOCT), [MaxBin2](https://sourceforge.net/projects/maxbin2/), and [MetaBAT2](https://sourceforge.net/projects/maxbin2/)
4. Refine & reassemble bins with [metaWRAP](https://github.com/bxlab/metaWRAP)
5. Taxonomic assignment with [GTDB-tk](https://github.com/Ecogenomics/GTDBTk)
6. Relative abundances with [bwa](https://github.com/lh3/bwa) and [samtools](https://github.com/samtools/samtools)
7. Reconstruct & evaluate genome-scale metabolic models with [CarveMe](https://github.com/cdanielmachado/carveme) and [memote](https://github.com/opencobra/memote)
8. Species metabolic coupling analysis with [SMETANA](https://github.com/cdanielmachado/smetana)

### Bonus

9. Growth rate estimation with [GRiD](https://github.com/ohlab/GRiD), [SMEG](https://github.com/ohlab/SMEG) or [CoPTR](https://github.com/tyjo/coptr)
10. Pangenome analysis with [roary](https://github.com/sanger-pathogens/Roary)
11. Eukaryotic draft bins with [EukRep](https://github.com/patrickwest/EukRep) and [EukCC](https://github.com/Finn-Lab/EukCC)

## 🏗️ Active Development

If you want to see any new additional or alternative tools incorporated into the `metaGEM` workflow please raise an issue or create a pull request. Snakemake allows workflows to be very flexible, so adding new rules is as easy as filling out the following template and adding it to the Snakefile:

```
rule package-name:
    input:
        rules.rulename.output
    output:
        f'{config["path"]["root"]}/{config["folder"]["X"]}/{{IDs}}/output.file'
    message:
        """
        Helpful and descriptive message detailing goal of this rule/package.
        """
    shell:
        """
        # Well documented command line instructions go here
        
        # Load conda environment 
        set +u;source activate {config[envs][package]};set -u;

        # Run tool
        package-name -i {input} -o {output}
        """
```

## 🖇️ Publications

The `metaGEM` workflow has been used in the following publications:

```
Plastic-degrading potential across the global microbiome correlates with recent pollution trends
J Zrimec, M Kokina, S Jonasson, F Zorrilla, A Zelezniak
MBio, 2021
```

```
Competition-cooperation in the chemoautotrophic ecosystem of Movile Cave: first metagenomic approach on sediments
Chiciudean, I., Russo, G., Bogdan, D.F. et al. 
Environmental Microbiome, 2022
```

```
The National Ecological Observatory Network’s soil metagenomes: assembly and basic analysis
Werbin ZR, Hackos B, Lopez-Nava J et al. 
F1000Research, 2022
```

[![arxiv](https://img.shields.io/badge/bioRxiv-10.1101%2F2020.12.13.422558%20-B31B1B)](https://www.biorxiv.org/content/10.1101/2020.12.13.422558v2.full)

## 🍾 Please cite

```
metaGEM: reconstruction of genome scale metabolic models directly from metagenomes
Francisco Zorrilla, Filip Buric, Kiran R Patil, Aleksej Zelezniak
Nucleic Acids Research, 2021; gkab815, https://doi.org/10.1093/nar/gkab815
``` 

[![Nucleic Acids Research](https://img.shields.io/badge/Nucleic%20Acids%20Research-10.1093%2Fnar%2Fgkab815-critical)](https://academic.oup.com/nar/advance-article/doi/10.1093/nar/gkab815/6382386)

## 📲 Contact

Please reach out with any comments, concerns, or discussions regarding `metaGEM`.

[![Gitter chat](https://badges.gitter.im/gitterHQ/gitter.png)](https://gitter.im/metaGEM/community)
[![Twitter](https://img.shields.io/badge/Twitter-%40metagenomez-lightblue)](https://twitter.com/metagenomez)
[![LinkedIn](https://img.shields.io/badge/LinkedIn-fzorrilla94-blue)](https://www.linkedin.com/in/fzorrilla94/)
[![email](https://img.shields.io/badge/email-fz274%40cam.ac.uk-%23a6bddb)](fz274@cam.ac.uk)
