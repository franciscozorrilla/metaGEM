[![Gitter chat](https://badges.gitter.im/gitterHQ/gitter.png)](https://gitter.im/metaGEM/community)

# `metaGEM` :gem:

> An easy-to-use workflow to generate context specific genome-scale metabolic models and predict metabolic interactions within microbial communities directly from metagenomic data.

![metawrapfigs_v2 002](https://user-images.githubusercontent.com/35606471/103545679-ceb71580-4e99-11eb-9862-084115121980.jpeg)

`metaGEM` integrates an array of existing bioinformatics and metabolic modeling tools using Snakemake, for the purpose of predicting metabolic interactions within bacterial communities of microbiomes. From whole metagenome shotgun datasets, metagenome assembled genomes (MAGs) are reconstructed, which are then converted into genome-scale metabolic models (GEMs) for *in silico* simulations of cross feeding interactions within sample based communities. Additional outputs include abundance estimates, taxonomic assignment, growth rate estimation, pangenome analysis, and eukaryotic MAG identification.

## :bulb: Quickstart

You can start using `metaGEM` with just one line of code:

```
git clone https://github.com/franciscozorrilla/metaGEM.git && cd metaGEM && rm -r .git && bash env_setup.sh
```

Congratulations, you can now start using `metaGEM`! Verify your installation by using the `check` task:

```
bash metaGEM.sh --task check
```

Please consult the [setup page](https://github.com/franciscozorrilla/metaGEM/wiki/Quickstart) in the wiki for further configuration instructions.

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

Snakefile wrapper/parser for metaGEM, for more details visit https://github.com/franciscozorrilla/metaGEM.

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

`metaGEM` can even be used to explore your own gut microbiome sequencing data from at-home-test-kit services such as [unseen bio](https://unseenbio.com/). The [following tutorial](https://github.com/franciscozorrilla/unseenbio_metaGEM) showcases the `metaGEM` workflow on two unseenbio samples.

## :book: Wiki

Refer to the [wiki](https://github.com/franciscozorrilla/metaGEM/wiki) for additional usage tips, frequently asked questions, and implementation details. 

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

## :heavy_check_mark: Please cite

```
metaGEM: reconstruction of genome scale metabolic models directly from metagenomes
Francisco Zorrilla, Kiran R. Patil, Aleksej Zelezniak
bioRxiv 2020.12.31.424982; doi: https://doi.org/10.1101/2020.12.31.424982 
```

## 📲 Contact

Please reach out with any comments, concerns, or discussions regarding `metaGEM`.

* Gitter: [@metaGEM](https://gitter.im/metaGEM/community)
* Twitter: [@metagenomez](https://twitter.com/metagenomez)
* LinkedIn: [@fzorrilla94](https://www.linkedin.com/in/fzorrilla94/)
* Email: fz274@cam.ac.uk

