#!/bin/bash

# Logo

printLogo() {

  echo '
                           /$$               /$$$$$$$   /$$$$$$   /$$$$$$            /$$                              
                          | $$              | $$__  $$ /$$__  $$ /$$__  $$          |__/                              
 /$$$$$$/$$$$   /$$$$$$  /$$$$$$    /$$$$$$ | $$  \ $$| $$  \ $$| $$  \__/  /$$$$$$  /$$  /$$$$$$   /$$$$$$   /$$$$$$$
| $$_  $$_  $$ /$$__  $$|_  $$_/   |____  $$| $$$$$$$ | $$$$$$$$| $$ /$$$$ /$$__  $$| $$ /$$__  $$ /$$__  $$ /$$_____/
| $$ \ $$ \ $$| $$$$$$$$  | $$      /$$$$$$$| $$__  $$| $$__  $$| $$|_  $$| $$  \ $$| $$| $$  \ $$| $$$$$$$$|  $$$$$$ 
| $$ | $$ | $$| $$_____/  | $$ /$$ /$$__  $$| $$  \ $$| $$  | $$| $$  \ $$| $$  | $$| $$| $$  | $$| $$_____/ \____  $$
| $$ | $$ | $$|  $$$$$$$  |  $$$$/|  $$$$$$$| $$$$$$$/| $$  | $$|  $$$$$$/| $$$$$$$/| $$| $$$$$$$/|  $$$$$$$ /$$$$$$$/
|__/ |__/ |__/ \_______/   \___/   \_______/|_______/ |__/  |__/ \______/ | $$____/ |__/| $$____/  \_______/|_______/ 
                                                                          | $$          | $$                          
                                                                          | $$          | $$                          
                                                                          |__/          |__/                          

A Snakemake-based metagenomics pipeline desinged to study microbial communities using high performance computer clusters.
'

}

# Helpfile function

usage() {

  printLogo

  echo -n "Usage: bash metaBAGpipes.sh [-t TASK] [-j NUMBER OF JOBS] [-c NUMBER OF CORES] [-m MEMORY]

Snakefile wrapper/parser for metaBAGpipes. 

 Options:
  -t, --task        Specify task to complete:

                        SETUP
                            createFolders
                            downloadToy
                            organizeData

                        WORKFLOW
                            fastp (Read QC)
                            megahit (Assembly)
                            kallisto (Mapping for CONCOCT)
                            concoct (Binning)
                            metabat (Binning)
                            maxbin (Binning)
                            binRefine (Bin refinement)
                            binReassemble (Bin reassembly)
                            classifyGenomes (Bin classification)
                            abundance (Bin quantification)
                            extractProteinBins
                            carveme
                            organizeGems
                            smetana
                            memote   
                            grid

                        VISUALIZATION (in development)
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

Example: bash metaBAGpipes.sh -t createFolders -j 1 -c 1

"
}

# Prepare to submit jobs function: unlock, dryrun, and display config files
snakePrep() {

        # Show config.yaml params
        echo -e "\nPlease verify parameters set in the config.yaml file: \n"
        paste config.yaml

        while true; do
            read -p "Do you wish to proceed with this config.yaml file? (y/n)" yn
            case $yn in
                [Yy]* ) echo " "; break;;
                [Nn]* ) exit;;
                * ) echo "Please answer yes or no.";;
            esac
        done

        # Show cluster_config.json params
        echo -e "Please verify parameters set in the cluster_config.json file: \n"
        paste cluster_config.json
        echo -e "\n"

        while true; do
            read -p "Do you wish to proceed with this cluster_config.json file? (y/n)" yn
            case $yn in
                [Yy]* ) echo " "; break;;
                [Nn]* ) exit;;
                * ) echo "Please answer yes or no.";;
            esac
        done

        echo -e "\nUnlocking snakemake ... "
        snakemake --unlock

        echo -e "\nDry-running snakemake jobs ... "
        snakemake all -j $njobs -n -k --cluster-config cluster_config.json -c "sbatch -A {cluster.account} -t {cluster.time} -n {cluster.n} --ntasks {cluster.tasks} --cpus-per-task {cluster.n} --output {cluster.output}"
}

# Submit login node function
submitLogin() {

    snakePrep

    while true; do
        read -p "Do you wish to submit this job? (y/n)" yn
        case $yn in
            [Yy]* ) snakemake $task ; break;;
            [Nn]* ) exit;;
            * ) echo "Please answer yes or no.";;
        esac
    done

}

# Submit cluster function

submitCluster() {

    # Parse Snakefile rule all (line 22 of Snakefile) input to match output of desired target rule stored in "$string". Note: Hardcoded line number.

    echo "Parsing Snakefile to target rule: $task ... "
    sed  -i "22s~^.*$~        $string~" Snakefile

    # Check if the number of cores flag is specified by user for cluster job

    if [[ -z "$ncores" ]]; then
        
        # No number of cores provided.
        echo -e "\nWARNING: User is requesting to submit cluster job without specifying the number of cores parameter (-n) ... "

    else

        # Parse cluster_config.json cores (line 5) to match number requested cores stored in "$ncores". Note: Hardcoded line number.

        echo "Parsing cluster_config.json to match requested number of cores: $ncores."
        sed -i "5s/:.*$/: $ncores,/" cluster_config.json

    fi 

    # Check if the number of jobs flag is specified by user for cluster job

    if [[ -z "$njobs" ]]; then
        
        # No number of jobs provided.
        echo "WARNING: User is requesting to submit cluster job without specifying the number of jobs parameter (-j) ... "

    fi   

    # Check if memory input argument was provided by user. If so, parse cluster_config.json memory (line 7) to match requested memory stored in "$mem". Note: Hardcoded line number.

    if [[ -z "$mem" ]]; then

        # No memory flag provided.
        echo "WARNING: User is requesting to submit cluster job without specifying the memory flag (-m) ... "

        snakePrep

        while true; do
            read -p "Do you wish to submit this batch of jobs? (y/n)" yn
            case $yn in
                [Yy]* ) echo "nohup snakemake all -j $njobs -k --cluster-config cluster_config.json -c 'sbatch -A {cluster.account} -t {cluster.time} -n {cluster.n} --ntasks {cluster.tasks} --cpus-per-task {cluster.n} --output {cluster.output}' &"|bash; break;;
                [Nn]* ) exit;;
                * ) echo "Please answer yes or no.";;
            esac
        done

    else

        # Memory flag was provided, parse cluster_config.json memory (line 7) to match number requested memory stored in "$mem". Note: Hardcoded line number.
        echo "Parsing cluster_config.json to match requested memory: $mem."
        sed -i "7s/:.*$/: $(echo $mem)G,/" cluster_config.json

        snakePrep

        while true; do
            read -p "Do you wish to submit this batch of jobs? (y/n)" yn
            case $yn in
                [Yy]* ) echo "nohup snakemake all -j $njobs -k --cluster-config cluster_config.json -c 'sbatch -A {cluster.account} -t {cluster.time} --mem {cluster.mem} -n {cluster.n} --ntasks {cluster.tasks} --cpus-per-task {cluster.n} --output {cluster.output}' &"|bash; break;;
                [Nn]* ) exit;;
                * ) echo "Please answer yes or no.";;
            esac
        done
    fi

}

# Parse function

parse() {

  printLogo

  # Set root folder

  echo "Setting current directory to root ... "
  root=$(pwd)
  sed  -i "2s~/.*$~$root~" config.yaml # hardcoded line for root, change the number 2 if any new lines are added to the start of config.yaml

  # No need to parse snakefile for login node jobs, submit the following locally

  if [ $task == "createFolders" ] || [ $task == "downloadToy" ] || [ $task == "organizeData" ] || [ $task == "qfilterVis" ] || [ $task == "assemblyVis" ] || [ $task == "binningVis" ] || [ $task == "taxonomyVis" ] ||  [ $task == "extractProteinBins" ] || [ $task == "organizeGems" ] || [ $task == "modelVis" ] || [ $task == "interactionVis" ] || [ $task == "growthVis" ]; then
    submitLogin

 # Parse snakefile for cluster jobs

  elif [ $task == "fastp" ]; then
    string='expand(config["path"]["root"]+"/"+config["folder"]["qfiltered"]+"/{IDs}/{IDs}_1.fastq.gz", IDs = IDs)'
    submitCluster

  elif [ $task == "megahit" ]; then
    string='expand(config["path"]["root"]+"/"+config["folder"]["assemblies"]+"/{IDs}/contigs.fasta.gz", IDs = IDs)'
    submitCluster

  elif [ $task == "kallisto" ]; then
    string='expand(config["path"]["root"]+"/"+config["folder"]["concoctInput"]+"/{IDs}_concoct_inputtableR.tsv", IDs = IDs)'
    submitCluster

  elif [ $task == "concoct" ]; then
    string='expand(config["path"]["root"]+"/"+config["folder"]["concoctOutput"]+"/{IDs}/{IDs}.concoct-bins", IDs = IDs)'
    submitCluster

  elif [ $task == "metabat" ]; then
    string='expand(config["path"]["root"]+"/"+config["folder"]["metabat"]+"/{IDs}/{IDs}.metabat-bins", IDs = IDs)'
    submitCluster

  elif [ $task == "maxbin" ]; then
    string='expand(config["path"]["root"]+"/"+config["folder"]["maxbin"]+"/{IDs}/{IDs}.maxbin-bins", IDs = IDs)'
    submitCluster

  elif [ $task == "binRefine" ]; then
    string='expand(config["path"]["root"]+"/"+config["folder"]["refined"]+"/{IDs}", IDs = IDs)'
    submitCluster

  elif [ $task == "binReassemble" ]; then
    string='expand(config["path"]["root"]+"/"+config["folder"]["reassembled"]+"/{IDs}", IDs = IDs)'
    submitCluster

  elif [ $task == "classifyGenomes" ]; then
    string='expand(config["path"]["root"]+"/"+config["folder"]["classification"]+"/{IDs}", IDs = IDs)'
    submitCluster

  elif [ $task == "abundance" ]; then
    string='expand(config["path"]["root"]+"/"+config["folder"]["abundance"]+"/{IDs}", IDs = IDs)'
    submitCluster

  elif [ $task == "carveme" ]; then
    string='expand(config["path"]["root"]+"/"+config["folder"]["GEMs"]+"/{binIDs}.xml", binIDs = binIDs)'
    submitCluster

  elif [ $task == "smetana" ]; then
    string='expand(config["path"]["root"]+"/"+config["folder"]["SMETANA"]+"/{IDs}_detailed.tsv", IDs = IDs)'
    submitCluster

  elif [ $task == "memote" ]; then
    string='expand(config["path"]["root"]+"/"+config["folder"]["memote"]+"/{IDs}", IDs = IDs)'
    submitCluster

  elif [ $task == "grid" ]; then
    string='expand(config["path"]["root"]+"/"+config["folder"]["GRiD"]+"/{IDs}", IDs = IDs)'
    submitCluster

  else
    echo "Task not recognized."
    usage
  fi

}

# Read input arguments

if [ $# -eq 0 ]; then
    echo "No arguments provided ... "
    usage
else
    # Read in options
    while [[ $1 = -?* ]]; do
      case $1 in
        -t|--task) shift; task=${1} ;;
        -j|--nJobs) shift; njobs=${1} ;;
        -c|--nCores) shift; ncores=${1} ;;
        -m|--mem) shift; mem=${1} ;;
        --endopts) shift; break ;;
        *) die "invalid option: '$1'." ;;
      esac
      shift
    done
    parse

fi
