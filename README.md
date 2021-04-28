# :gem: `metaGEM`

> An easy-to-use workflow for generating context specific genome-scale metabolic models and predicting metabolic interactions within microbial communities directly from metagenomic data.

[![DOI](https://img.shields.io/badge/bioRxiv-10.1101%2F2020.12.31.424982%20-B31B1B)](https://www.biorxiv.org/content/10.1101/2020.12.31.424982v2.full)
[![Build Status](https://travis-ci.org/franciscozorrilla/metaGEM.svg?branch=master)](https://travis-ci.org/franciscozorrilla/metaGEM)
[![GitHub license](https://img.shields.io/github/license/franciscozorrilla/metaGEM)](https://github.com/franciscozorrilla/metaGEM/blob/master/LICENSE)
[![Gitter chat](https://badges.gitter.im/gitterHQ/gitter.png)](https://gitter.im/metaGEM/community)
[![DOI](https://img.shields.io/badge/Zenodo-10.5281%2F4707723-blue)](https://zenodo.org/badge/latestdoi/137376259)
[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/drive/1I1S8AoGuJ9Oc2292vqAGTDmZcbnolbuj#scrollTo=awiAaVwSF5Fz)

![metawrapfigs_final 001](https://user-images.githubusercontent.com/35606471/114914819-b2246500-9e1a-11eb-8daa-3860ec92af30.jpeg)

`metaGEM` integrates an array of existing bioinformatics and metabolic modeling tools using Snakemake, for the purpose of predicting metabolic interactions within bacterial communities of microbiomes. From whole metagenome shotgun datasets, metagenome assembled genomes (MAGs) are reconstructed, which are then converted into genome-scale metabolic models (GEMs) for *in silico* simulations. Additional outputs include abundance estimates, taxonomic assignment, growth rate estimation, pangenome analysis, and eukaryotic MAG identification. The workflow diagram highlights the most important tasks in the core workflow on the left-hand column, whiel the auxiliary or supplementary tools are displayed on the right-hand column.

## 🌥️ Try it now

You can analyze a toy dataset using `metaGEM` by following along the google colab notebook on the cloud.

[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/drive/1I1S8AoGuJ9Oc2292vqAGTDmZcbnolbuj#scrollTo=awiAaVwSF5Fz)

## :bulb: Installation

You can set up `metaGEM` on your cluster with just one line of code:

```
git clone https://github.com/franciscozorrilla/metaGEM.git && cd metaGEM && rm -r .git && bash env_setup.sh
```

Congratulations, you can now start using `metaGEM`! Verify your installation by using the `check` task:

```
bash metaGEM.sh --task check
```

Please consult the setup page in the wiki for further configuration instructions.

[![installation](https://img.shields.io/badge/metaGEM-Installation-%2331a354)](https://github.com/franciscozorrilla/metaGEM/wiki/Quickstart)

## 🔧 Usage

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

## 💩 Tutorial

`metaGEM` can even be used to explore your own gut microbiome sequencing data from at-home-test-kit services such as [unseen bio](https://unseenbio.com/). The following tutorial showcases the `metaGEM` workflow on two unseenbio samples.

[![Tutorial](https://img.shields.io/badge/metaGEM-Tutorial-%23d8b365)](https://github.com/franciscozorrilla/unseenbio_metaGEM)

## :book: Wiki

Refer to the wiki for additional usage tips, frequently asked questions, and implementation details.

[![wiki](https://img.shields.io/badge/metaGEM-Wiki-blue)](https://github.com/franciscozorrilla/metaGEM/wiki)

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

## :construction: Active Development

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

## :paperclip: Publications

The `metaGEM` workflow was used in the following publication(s):

```
Plastic-degrading potential across the global microbiome correlates with recent pollution trends
Jan Zrimec, Mariia Kokina, Sara Jonasson, Francisco Zorrilla, Aleksej Zelezniak
bioRxiv 2020.12.13.422558; doi: https://doi.org/10.1101/2020.12.13.422558 
```

[![arxiv](https://img.shields.io/badge/bioRxiv-10.1101%2F2020.12.13.422558%20-B31B1B)](https://www.biorxiv.org/content/10.1101/2020.12.13.422558v2.full)

## :heavy_check_mark: Please cite

```
metaGEM: reconstruction of genome scale metabolic models directly from metagenomes
Francisco Zorrilla, Kiran R. Patil, Aleksej Zelezniak
bioRxiv 2020.12.31.424982; doi: https://doi.org/10.1101/2020.12.31.424982 
```
[![DOI](https://img.shields.io/badge/bioRxiv-10.1101%2F2020.12.31.424982%20-B31B1B)](https://www.biorxiv.org/content/10.1101/2020.12.31.424982v2.full)

## 📲 Contact

Please reach out with any comments, concerns, or discussions regarding `metaGEM`.

[![Gitter chat](https://badges.gitter.im/gitterHQ/gitter.png)](https://gitter.im/metaGEM/community)
[![Twitter](https://img.shields.io/badge/Twitter-%40metagenomez-lightblue)](https://twitter.com/metagenomez)
[![LinkedIn](https://img.shields.io/badge/LinkedIn-fzorrilla94-blue)](https://www.linkedin.com/in/fzorrilla94/)
[![email](https://img.shields.io/badge/email-fz274%40cam.ac.uk-%23a6bddb)](fz274@cam.ac.uk)
