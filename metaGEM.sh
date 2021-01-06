#!/bin/bash

# Logo
printLogo() {

  echo '

_________________________________________________________________________/\\\\\\\\\\\\___/\\\\\\\\\\\\\\\___/\\\\____________/\\\\_        
 _______________________________________________________________________/\\\//////////___\/\\\///////////___\/\\\\\\________/\\\\\\_       
  __________________________________________/\\\________________________/\\\______________\/\\\______________\/\\\//\\\____/\\\//\\\_      
   ____/\\\\\__/\\\\\________/\\\\\\\\____/\\\\\\\\\\\___/\\\\\\\\\_____\/\\\____/\\\\\\\__\/\\\\\\\\\\\______\/\\\\///\\\/\\\/_\/\\\_     
    __/\\\///\\\\\///\\\____/\\\/////\\\__\////\\\////___\////////\\\____\/\\\___\/////\\\__\/\\\///////_______\/\\\__\///\\\/___\/\\\_    
     _\/\\\_\//\\\__\/\\\___/\\\\\\\\\\\______\/\\\_________/\\\\\\\\\\___\/\\\_______\/\\\__\/\\\______________\/\\\____\///_____\/\\\_   
      _\/\\\__\/\\\__\/\\\__\//\\///////_______\/\\\_/\\____/\\\/////\\\___\/\\\_______\/\\\__\/\\\______________\/\\\_____________\/\\\_  
       _\/\\\__\/\\\__\/\\\___\//\\\\\\\\\\_____\//\\\\\____\//\\\\\\\\/\\__\//\\\\\\\\\\\\/___\/\\\\\\\\\\\\\\\__\/\\\_____________\/\\\_ 
        _\///___\///___\///_____\//////////_______\/////______\////////\//____\////////////_____\///////////////___\///______________\///__


A Snakemake-based metagenomics pipeline desinged to study the metabolism of microbial communities using high performance computer clusters.
'

}

# Helpfile function

usage() {

  printLogo

  echo -n "Usage: bash metaGEM.sh [-t|--task TASK] 
                [-j|--nJobs NUMBER OF JOBS] 
                [-c|--cores NUMBER OF CORES] 
                [-m|--mem GB RAM] 
                [-h|--hours MAX RUNTIME]

Snakefile wrapper/parser for metaGEM. 

 Options:
  -t, --task        Specify task to complete:

                        SETUP
                            createFolders
                            downloadToy
                            organizeData

                        WORKFLOW
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
                            grid
                            prokka
                            roary

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
  -h, --hours       Specify number of hours to allocated to job runtime

Suggested workflow:

    0. metaGEM setup
    1. Quality filter reads with fastp
    2. Assembly with megahit
    3. Draft bin sets with CONCOCT,MaxBin2, and MetaBAT2
    4. Refine & reassemble bins with metaWRAP
    5. Taxonomic assignment with GTDB-tk
    6. Relative abundances with bwa and samtools
    7. Reconstruct & evaluate genome-scale metabolic models with CarveMe and memote
    8. Species metabolic coupling analysis with SMETANA
    9. Growth rate estimation with GRiD
    10. Pangenome analysis with roary
    11. Eukaryotic draft bins with EukRep and EukCC


e.g. to submit 10 short read quality filtering jobs with 2 cores + 4 GB RAM each and maximum runtime of 1 hour:
     bash metaGEM.sh -t fastp -j 10 -c 2 -m 4 -h 1

"
}

# Prompt user to confirm input parameters/options
checkParams() {

    echo " "
    while true; do
        read -p "Do you wish to continue with these parameters? (y/n)" yn
        case $yn in
            [Yy]* ) echo "Proceeding with $task job(s) ... " ; break;;
            [Nn]* ) exit;;
            * ) echo "Please answer yes or no.";;
        esac
    done

}

# Display config.yaml function for user inspection
snakeConfig() {

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

}

# Display cluster_config.json function for user inspection
clusterConfig() {

    # Show cluster_config.json params
    echo -e "Please verify parameters set in the cluster_config.json file: \n"
    paste cluster_config.json
    echo " "

    while true; do
        read -p "Do you wish to proceed with this cluster_config.json file? (y/n)" yn
        case $yn in
            [Yy]* ) echo " "; break;;
            [Nn]* ) exit;;
            * ) echo "Please answer yes or no.";;
        esac
    done
    
}

# Prepare to submit cluster jobs function: display config files, unlock, and dry run
snakePrep() {

    snakeConfig

    clusterConfig

    echo "Unlocking snakemake ... "
    snakemake --unlock -j 1

    echo -e "\nDry-running snakemake jobs ... "
    snakemake all -j $njobs -n -k --cluster-config cluster_config.json -c "sbatch -A {cluster.account} -t {cluster.time} -n {cluster.n} --ntasks {cluster.tasks} --cpus-per-task {cluster.n} --output {cluster.output}"
}

# Submit login node function
submitLogin() {

    echo "No need to parse Snakefile for target rule: $task ... "

    checkParams

    snakeConfig

    echo "Unlocking snakemake ... "
    snakemake --unlock -j 1
    echo " "

    while true; do
        read -p "Do you wish to submit this $task job? (y/n)" yn
        case $yn in
            [Yy]* ) snakemake $task -j 1; break;;
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

    # Check if the number of jobs flag is specified by user for cluster job
    if [[ -z "$njobs" ]]; then
        
        # No number of jobs provided.
        echo "WARNING: User is requesting to submit cluster job without specifying the number of jobs parameter (-j) ... "
    else   
        # Number of jobs provided.
        echo "Number of jobs to be sumitted to cluster: $njobs ... "
    fi   

    # Check if the number of cores flag is specified by user for cluster job
    if [[ -z "$ncores" ]]; then
        
        # No number of cores provided.
        echo "WARNING: User is requesting to submit cluster job without specifying the number of cores parameter (-c) ... "

    else

        # Parse cluster_config.json cores (line 5) to match number requested cores stored in "$ncores". Note: Hardcoded line number.
        echo "Parsing cluster_config.json to match requested number of cores: $ncores ... "
        sed -i "5s/:.*$/: $ncores,/" cluster_config.json

    fi 

    # Check if the hours flag is specified by user for cluster job
    if [[ -z "$hours" ]]; then
        
        # No number of jobs provided.
        echo "WARNING: User is requesting to submit cluster job without specifying the number of hours parameter (-h) ... "

    else

        # Parse cluster_config.json time (line 4) to match number requested hours stored in "$hours". Note: Hardcoded line number.
        echo "Parsing cluster_config.json to match requested time (hours): $hours ... "
        sed -i "4s/:.*$/: \"0-$hours:00:00\",/" cluster_config.json

    fi 

    # Check if memory input argument was provided by user. If so, parse cluster_config.json memory (line 7) to match requested memory stored in "$mem". Note: Hardcoded line number.
    if [[ -z "$mem" ]]; then

        # No memory flag provided.
        echo "WARNING: User is requesting to submit cluster job without specifying the memory flag (-m) ... "

        checkParams

        snakePrep

        while true; do
            read -p "Do you wish to submit this batch of $task jobs? (y/n)" yn
            case $yn in
                [Yy]* ) echo "nohup snakemake all -j $njobs -k --cluster-config cluster_config.json -c 'sbatch -A {cluster.account} -t {cluster.time} -n {cluster.n} --ntasks {cluster.tasks} --cpus-per-task {cluster.n} --output {cluster.output}' &"|bash; break;;
                [Nn]* ) exit;;
                * ) echo "Please answer yes or no.";;
            esac
        done

    else

        # Memory flag was provided, parse cluster_config.json memory (line 7) to match number requested memory stored in "$mem". Note: Hardcoded line number.
        echo "Parsing cluster_config.json to match requested memory: $mem ... "
        sed -i "7s/:.*$/: $(echo $mem)G,/" cluster_config.json

        checkParams

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
  if [ $task == "createFolders" ] || [ $task == "downloadToy" ] || [ $task == "organizeData" ] || [ $task == "qfilterVis" ] || [ $task == "assemblyVis" ] || [ $task == "binningVis" ] || [ $task == "taxonomyVis" ] || [ $task == "abundanceVis" ] || [ $task == "extractProteinBins" ] || [ $task == "extractDnaBins" ] || [ $task == "organizeGEMs" ] || [ $task == "modelVis" ] || [ $task == "interactionVis" ] || [ $task == "growthVis" ] || [ $task == "prepareRoary" ]; then
    submitLogin

 # Parse snakefile for cluster jobs
  elif [ $task == "fastp" ]; then
    string='expand(config["path"]["root"]+"/"+config["folder"]["qfiltered"]+"/{IDs}/{IDs}_1.fastq.gz", IDs = IDs)'
    submitCluster

  elif [ $task == "megahit" ]; then
    string='expand(config["path"]["root"]+"/"+config["folder"]["assemblies"]+"/{IDs}/contigs.fasta.gz", IDs = IDs)'
    submitCluster

  elif [ $task == "crossMap" ]; then
    string='expand(config["path"]["root"]+"/"+config["folder"]["concoct"]+"/{IDs}/cov", IDs = IDs)'
    submitCluster

  elif [ $task == "concoct" ]; then
    string='expand(config["path"]["root"]+"/"+config["folder"]["concoct"]+"/{IDs}/{IDs}.concoct-bins", IDs = IDs)'
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

  elif [ $task == "gtdbtk" ]; then
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
    string='expand(config["path"]["root"]+"/"+config["folder"]["memote"]+"/{gemIDs}", gemIDs = gemIDs)'
    submitCluster

  elif [ $task == "grid" ]; then
    string='expand(config["path"]["root"]+"/"+config["folder"]["GRiD"]+"/{IDs}", IDs = IDs)'
    submitCluster

  elif [ $task == "prokka" ]; then
    string='expand(config["path"]["root"]+"/"+config["folder"]["pangenome"]+"/prokka/unorganized/{binIDs}", binIDs = binIDs)'
    submitCluster

  elif [ $task == "roary" ]; then
    string='expand(config["path"]["root"]+"/"+config["folder"]["pangenome"]+"/roary/{speciesIDs}/", speciesIDs = speciesIDs)'
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
        -h|--hours) shift; hours=${1} ;;
        --endopts) shift; break ;;
        * ) echo "Unknown option(s) provided, please read helpfile ... " && usage && exit 1;;
      esac
      shift
    done
    parse

fi
