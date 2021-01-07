#!/bin/bash

echo -e "\n==============================================\n
    		metaGEM setup\n
==============================================\n"

echo -e "Creating metagem conda environment and installing tools ... \n"

conda env create -f metaGEM_env.yml

echo -e "Downloading GTDB-tk database, requires ~25 Gb ... \n"

download-db.sh

echo -e "Creating metawrap conda environment and installing tools ... \n"

conda env create -f metaWRAP_env.yml

echo -e "Creating prokka + roary conda environment and installing tools ... \n"

conda env create -f prokkaroary_env.yml

echo -e "Please ensure that the installation directory is present in your $PATH variable if any installation issues arise."
