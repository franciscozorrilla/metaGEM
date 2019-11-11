## metaBAGpipes: metaBolism And metaGenomic analysis pipelines

### Use cases

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


### Installation

```bash
git clone git@github.com:franciscozorrilla/metaBAGpipes.git
```

#### Conda

A [conda](https://conda.io/en/latest/) specification file is provided to install required packages inside an
environment called `metabagpipes`.

```bash
conda env create -f metaBAGpipes_env.yml
source activate metabagpipes
```

#### Singularity

A [Singularity](https://sylabs.io/docs/) recipe files is provided to build an image that can be used with HPC clusters.

```bash
sudo singularity --verbose build metabagpipes.simg Singularity
```

### Toy dataset

As a tutorial, and to verify that your metaBAGpipes installation is working correctly, we provide a [toy dataset](https://zenodo.org/record/3534949#.XclQriV7lTZ). Expected output can be found in the Tutorial folder.

### Abstract
metaBAGpipes integrates an array of existing bioinformatics and metabolic modeling tools using Snakemake, for the purpose of interrogating social interactions in bacterial communities of the human gut microbiome. From WGS metagenomic datasets, metagenome assembled genomes (MAGs) are reconstructed, which are then converted into genome-scale metabolic models (GEMs) for *in silico* simulations of cross feeding interactions within sample based communities. Abundance estimates for community members are estimated by mapping metagenomic samples to the generated MAGs, which are used in combination with the simulated cross feeding interactions for the generation of explanatory and statistically significant linear models. We conclude that there is indeed a correlation, ranging from weak to moderate, between gut microbiome members’ abundance and set of metabolic cross-feeding interactions across samples. A more comprehensive analysis incorporating multiple datasets needs to be conducted to strengthen and expand the findings of this work.

# ![pipemap_v0.1.png](pipemap_v0.1.png)

### Significance

To our best knowledge, the work presented in this thesis project represents the first application of sample-specific gut microbiome community FBA simulations using MAG-based GEMs. While there are limitations to this approach, as discussed in the thesis write-up, it is conceivable that a fine-tuned version of such a pipeline could be used to study, develop, and evualuate personalized medicine treatments. 137 WGS gut metagenome samples were processed in this study, a relatively modest amount of data. By analyzing a greater number of samples, further machine learning methods can be employed to interrogate the network architechture of microbial communities and facilitiate data driven hypothesis generation. The figure below is a visual representation of the output information one can obtain from a particular sample using metaBAGpipes.

# ![ERR260149.png](ERR260149.png)

Dataset used:
  * Karlsson, Fredrik H., et al. “Gut Metagenome in European Women with Normal, Impaired and Diabetic Glucose Control.” *Nature*, vol.498, no.7452, 2013, pp.99–103. , doi:10.1038/nature12198.

This repository is administered by Francisco Zorrilla ([@franciscozorrilla](https://github.com/franciscozorrilla/)), Structural and Computational Biology Unit, EMBL. metaBAGpipes was developed throughout my Master's thesis project at the Systems and Synthetic Biology division of Chalmers Univeristy of Technology, under the supervision of Aleksej Zelezniak.

  * Last update: 11-11-2019
