#!/bin/bash

echo -e "\n==============================================================================================================\n
        		    			metaGEM setup\n
==============================================================================================================\n"

while true; do
    read -p "Do you wish to download and set up metaGEM conda environment? (y/n)" yn
    case $yn in
        [Yy]* ) echo "conda env create -f metaGEM_env.yml && source activate metagem && pip install --user memote carveme smetana && conda deactivate"|bash; break;;
        [Nn]* ) echo -e "\nSkipping metaGEM env setup, note that you will need this for refinement & reassembly of MAGs.\n"; break;;
        * ) echo "Please answer yes or no.";;
    esac
done

while true; do
    read -p "Do you wish to download the GTDB-tk database (~25 Gb)? (y/n)" yn
    case $yn in
        [Yy]* ) echo "download-db.sh"|bash; break;;
        [Nn]* ) echo -e "\nSkipping GTDB-tk database download, note that you will need this for taxonomic classification of MAGs.\n"; break;;
        * ) echo "Please answer yes or no.";;
    esac
done

while true; do
    read -p "Do you wish to download and set up metaWRAP conda environment? (y/n)" yn
    case $yn in
        [Yy]* ) echo "conda env create -f metaWRAP_env.yml"|bash; break;;
        [Nn]* ) echo -e "\nSkipping metaWRAP env setup, note that you will need this for refinement & reassembly of MAGs.\n"; break;;
        * ) echo "Please answer yes or no.";;
    esac
done

while true; do
    read -p "Do you wish to download the CheckM database (~275 Mb)? (y/n)" yn
    case $yn in
        [Yy]* ) echo "wget https://data.ace.uq.edu.au/public/CheckM_databases/checkm_data_2015_01_16.tar.gz"|bash; break;;
        [Nn]* ) echo -e "\nSkipping CheckM database download, note that you will need this for bin refinement & reassembly.\n"; break;;
        * ) echo "Please answer yes or no.";;
    esac
done

while true; do
    read -p "Do you wish to download and set up prokka + roary conda environment? (y/n)" yn
    case $yn in
        [Yy]* ) echo "conda env create -f prokkaroary_env.yml"|bash; break;;
        [Nn]* ) echo -e "\nSkipping CheckM database download, note that you will need this for bin refinement & reassembly.\n"; break;;
        * ) echo "Please answer yes or no.";;
    esac
done

echo 'Please ensure that the installation directory is present in your $PATH variable if installation issues arise with any tools.'
echo ""