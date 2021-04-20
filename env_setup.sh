#!/bin/bash

  echo '
=================================================================================================================================
Developed by: Francisco Zorrilla, Kiran R. Patil, and Aleksej Zelezniak___________________________________________________________
Publication: doi.org/10.1101/2020.12.31.424982___________________________/\\\\\\\\\\\\___/\\\\\\\\\\\\\\\___/\\\\____________/\\\\_        
________________________________________________________________________/\\\//////////___\/\\\///////////___\/\\\\\\________/\\\\\\_       
____________________________________________/\\\________________________/\\\______________\/\\\______________\/\\\//\\\____/\\\//\\\_      
_______/\\\\\__/\\\\\________/\\\\\\\\____/\\\\\\\\\\\___/\\\\\\\\\_____\/\\\____/\\\\\\\__\/\\\\\\\\\\\______\/\\\\///\\\/\\\/_\/\\\_     
______/\\\///\\\\\///\\\____/\\\/////\\\__\////\\\////___\////////\\\____\/\\\___\/////\\\__\/\\\///////_______\/\\\__\///\\\/___\/\\\_    
______\/\\\_\//\\\__\/\\\___/\\\\\\\\\\\______\/\\\_________/\\\\\\\\\\___\/\\\_______\/\\\__\/\\\______________\/\\\____\///_____\/\\\_   
_______\/\\\__\/\\\__\/\\\__\//\\///////_______\/\\\_/\\____/\\\/////\\\___\/\\\_______\/\\\__\/\\\______________\/\\\_____________\/\\\_  
________\/\\\__\/\\\__\/\\\___\//\\\\\\\\\\_____\//\\\\\____\//\\\\\\\\/\\__\//\\\\\\\\\\\\/___\/\\\\\\\\\\\\\\\__\/\\\_____________\/\\\_ 
_________\///___\///___\///_____\//////////_______\/////______\////////\//____\////////////_____\///////////////___\///______________\///__
=============================================================================================================================================
A Snakemake-based pipeline desinged to predict metabolic interactions directly from metagenomics data using high performance computer clusters
===============================================================================================================================================
'

#check if conda is installed/available
echo -ne "Checking if conda is available ... "
condatest=$(conda list|wc -l)

if [[ "$condatest" -eq 0 ]]; then
    echo -e "\nWARNING: Conda is not available! Please load your cluster's conda module or install locally and re-run the env_setup.sh script using:\n\nbash env_setup.sh\n" && exit
elif [[ "$condatest" -gt 0 ]]; then
    condav=$(conda --version|cut -d ' ' -f2)
    echo -e "detected version $condav!"
fi

while true; do
    read -p "Do you wish to download and set up metaGEM conda environment? (y/n)" yn
    case $yn in
        [Yy]* ) echo "conda create -n mamba mamba && source activate mamba && mamba env create -f envs/metaGEM_env.yml && pip install --user memote carveme smetana && conda deactivate && echo "|bash; break;;
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
        [Yy]* ) echo "conda env create -f envs/metaWRAP_env.yml"|bash; break;;
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
        [Yy]* ) echo "conda env create -f envs/prokkaroary_env.yml"|bash; break;;
        [Nn]* ) echo -e "\nSkipping prokka-roary env setup, note that you will need this for pangenome analysis of MAGs.\n"; break;;
        * ) echo "Please answer yes or no.";;
    esac
done

echo 'Please ensure that the installation directory is present in your $PATH variable if installation issues arise with any tools.'
echo ""