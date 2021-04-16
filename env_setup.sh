#!/bin/bash

echo "========================================================================================================================================
        		    			
_________________________________________________________________________/\\\\\\\\\\\\___/\\\\\\\\\\\\\\\___/\\\\____________/\\\\_        
 _______________________________________________________________________/\\\//////////___\/\\\///////////___\/\\\\\\________/\\\\\\_       
  __________________________________________/\\\________________________/\\\______________\/\\\______________\/\\\//\\\____/\\\//\\\_      
   ____/\\\\\__/\\\\\________/\\\\\\\\____/\\\\\\\\\\\___/\\\\\\\\\_____\/\\\____/\\\\\\\__\/\\\\\\\\\\\______\/\\\\///\\\/\\\/_\/\\\_     
    __/\\\///\\\\\///\\\____/\\\/////\\\__\////\\\////___\////////\\\____\/\\\___\/////\\\__\/\\\///////_______\/\\\__\///\\\/___\/\\\_    
     _\/\\\_\//\\\__\/\\\___/\\\\\\\\\\\______\/\\\_________/\\\\\\\\\\___\/\\\_______\/\\\__\/\\\______________\/\\\____\///_____\/\\\_   
      _\/\\\__\/\\\__\/\\\__\//\\///////_______\/\\\_/\\____/\\\/////\\\___\/\\\_______\/\\\__\/\\\______________\/\\\_____________\/\\\_  
       _\/\\\__\/\\\__\/\\\___\//\\\\\\\\\\_____\//\\\\\____\//\\\\\\\\/\\__\//\\\\\\\\\\\\/___\/\\\\\\\\\\\\\\\__\/\\\_____________\/\\\_ 
        _\///___\///___\///_____\//////////_______\/////______\////////\//____\////////////_____\///////////////___\///______________\///__

========================================================================================================================================"

while true; do
    read -p "Do you wish to download and set up metaGEM conda environment? (y/n)" yn
    case $yn in
        [Yy]* ) echo "conda env create -f envs/metaGEM_env.yml && source activate metagem && pip install --user memote carveme smetana && conda deactivate"|bash; break;;
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