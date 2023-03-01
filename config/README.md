# 💎 Setup guide

## 🔩 Config files

Make sure to inspect and set up the two config files in this folder.

### Snakemake configuration 
`config.yaml`: handles all the tunable parameters, subfolder names, paths, and more. The `root` path is automatically set by the `metaGEM.sh` parser to be the current working directory. Most importantly, you should make sure that the `scratch` path is properly configured. Most clusters have a location for temporary or high I/O operations such as `$TMPDIR` or `$SCRATCH`, e.g. [see here](https://www.c3se.chalmers.se/documentation/filesystem/#using-node-local-disk-tmpdir). Please refer to the config.yaml [wiki page](https://github.com/franciscozorrilla/metaGEM/wiki/Snakefile-config) for a more in depth look at this config file.

### Cluster configuration
`cluster_config.json`: handles parameters for submitting jobs to the cluster workload manager. Most importantly, you should make sure that the `account` is properly defined to be able to submit jobs to your cluster. Please refer to the cluster_config.json wiki page for a more in depth look at this config file.

## 🛢️ Environments

Set up three conda environments:
1. `mamba`: Used for installing mamba and setting up subsequent environments from recipe files
2. `metagem`: Contains most `metaGEM` core workflow tools, Python 3
3. `metawrap` Contains `metaWRAP` and its dependencies, Python 2

### 1. mamba

Conda can take *ages* to solve environment dependencies when installing many tools at once, we can use [mamba](https://github.com/mamba-org/mamba) instead for faster installation.

```
conda create -n mamba mamba
```

Activate mamba environment to quickly set up subsequent environments.

```
source activate mamba
```

### 2. metaGEM

Clone metaGEM repo

```
git clone https://github.com/franciscozorrilla/metaGEM.git
```

Move into metaGEM/workflow folder

```
cd metaGEM/workflow
```

Clean up unnecessary ~250 Mb of unnecessary git history files

```
rm -r ../.git
```

Press `y` and `Enter` when prompted to remove write-protected files, these are not necessary and just eat your precious space.

```
rm: remove write-protected regular file ‘.git/objects/pack/pack-f4a65f7b63c09419a9b30e64b0e4405c524a5b35.pack’? y
rm: remove write-protected regular file ‘.git/objects/pack/pack-f4a65f7b63c09419a9b30e64b0e4405c524a5b35.idx’? y
```

Recommended method: install from bioconda
```
conda config --add channels conda-forge && mamba create --prefix envs/metagem -c bioconda metagem 
```

Alternative method: create metaGEM env using recipe .yml file

```
mamba env create --prefix ./envs/metagem -f envs/metaGEM_env.yml
```

Deactivate mamba env and activate metaGEM env

```
source deactivate && source activate envs/metagem
```

Install pip tools

```
pip install --user memote carveme smetana
```

### 3. metaWRAP

It is best to set up `metaWRAP` in its own isolated environment to prevent version conflicts with `metaGEM`. Note that `metaWRAP v1.3.2` has not migrated from python 2 to python 3 yet.

```
conda create -n metawrap
source activate metawrap
conda install -c ursky metawrap-mg=1.3.2
```

Or using the conda recipe file:

```
mamba env create --prefix ./envs/metawrap -f envs/metaWRAP_env.yml
```

## 🔮 Check installation

To make sure that the basics have been properly configured, run the `check` task using the `metaGEM.sh` parser:

```
bash metaGEM.sh -t check
```

This will check if conda is installed/available and verify that the environments were properly set up.
Additionally, this `check` function will prompt you to create results folders if they are not already present.
Finally, this task will check if any sequencing files are present in the dataset folder, prompting the user to the either organize already existing files into sample-specific subfolders or to download a small [toy dataset](https://zenodo.org/record/3534949/). 

## Tools requiring additional configuration

> **Warning** Please note that you will need to set up the following tools/databases to run the complete core metaGEM workflow:

### 1. CheckM

`CheckM` is used extensively within the `metaWRAP` modules to evaluate the output of various intermediate steps. Although the `CheckM` package is installed in the `metawrap` environment, the user is required to download the `CheckM` database and run `checkm data setRoot <db_dir>` as outlined in the [`CheckM` installation guide](https://github.com/Ecogenomics/CheckM/wiki/Installation#how-to-install-checkm).

### 2. GTDB-Tk

`GTDB-Tk` is used for taxonomic assignment of MAGs, and requires a database to be downloaded and configured. Please refer to the [installation documentation](https://ecogenomics.github.io/GTDBTk/installing/index.html) for detailed instructions.

### 3. CPLEX

Unfortunately `CPLEX` cannot be automatically installed in the `env_setup.sh` script, you must install this dependency manually within the metagem conda environment. GEM reconstruction and GEM community simulations require the `IBM CPLEX solver`, which is [free to download with an academic license](https://www.ibm.com/academic/home). Refer to the [`CarveMe`](https://carveme.readthedocs.io/en/latest/installation.html) and [`SMETANA`](https://smetana.readthedocs.io/en/latest/installation.html) installation instructions for further information or troubleshooting. Note: `CPLEX v.12.8` is recommended.
