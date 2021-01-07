#!/bin/bash

echo -e "\n==============================================\n
    		metaGEM setup\n
==============================================\n"

echo -e "Creating metagem conda environment and installing tools ... \n"
conda env create -f metaGEM_env.yml

echo -e "Installing remaining tools via pip ... \n"
source activate metagem
pip install --user memote carveme smetana

echo -e "Downloading GTDB-tk database, requires ~25 Gb ... \n"
download-db.sh

echo -e "Done creating metagem conda environment!\n"
conda deactivate

echo -e "Creating metawrap conda environment and installing tools ... \n"
conda env create -f metaWRAP_env.yml

echo -e "Downloading CheckM database, requires ~275 Mb ... \n"
wget https://data.ace.uq.edu.au/public/CheckM_databases/checkm_data_2015_01_16.tar.gz

echo -e "Done creating metawrap conda environment!\n"

echo -e "Creating prokka + roary conda environment and installing tools ... \n"
conda env create -f prokkaroary_env.yml

echo -e "Done creating prokka + roary conda environment!\n"

echo -e "Please ensure that the installation directory is present in your $PATH variable if any installation issues arise."