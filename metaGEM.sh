#!/bin/bash

# Logo
printLogo() {

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

}

# Helpfile function

usage() {

  printLogo

  echo -n "Usage: bash metaGEM.sh [-t|--task TASK] 
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
                            crossMapSeries
                            kallistoIndex
                            crossMapParallel
                            kallisto2concoct
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

"
}

# Run check task
run_check() {

#check if conda is installed/available
echo -ne "Checking if conda is available ... "
condatest=$(conda list|wc -l)
if [[ "$condatest" -eq 0 ]]; then
    echo -e "WARNING: Conda is not available! Please load your cluster's conda module or install locally.\n" && exit
elif [[ "$condatest" -gt 0 ]]; then
    condav=$(conda --version|cut -d ' ' -f2)
    echo -e "detected version $condav!"
fi

# check if conda environments are present
echo -ne "Searching for metaGEM conda environment ... "
envcheck1=$(conda info --envs|grep -w metagem|wc -l)
if [[ "$envcheck1" -ge 1 ]]; then
    echo "detected! Activating metagem env ... "
    source activate metagem
else
    echo "not detected, please run the env_setup.sh script!"
fi

echo -ne "Searching for metaWRAP conda environment ... "
envcheck2=$(conda info --envs|grep -w metawrap|wc -l)
if [[ "$envcheck2" -ge 1 ]]; then
    echo "detected!"
else
    echo "not detected, please run the env_setup.sh script!"
fi

echo -ne "Searching for prokka-roary conda environment ... "
envcheck3=$(conda info --envs|grep -w prokkaroary|wc -l)
if [[ "$envcheck3" -ge 1 ]]; then
    echo -e "detected!\n"
else
    echo -e "not detected, please run the env_setup.sh script!\n"
fi

# run createFolders rule to create folders in case any of them are missing
echo -e "Checking folders in workspace $pwd ... "
nFolders=$(ls -d */|wc -l)
if [[ "$nFolders" -le 20 ]]; then
    while true; do
        read -p "Some folders appear to be missing, do you wish to run the createFolders Snakefile rule? (y/n)" yn
        case $yn in
            [Yy]* ) echo "Running the createFolders snakefile rule ... " && snakemake createFolders -j1; break;;
            [Nn]* ) echo "Skipping folder creation ... "; break;;
            * ) echo "Please answer yes or no.";;
        esac
    done
fi

# search for folders and files with .gz extension within dataset folder
count_files=$(find dataset -name "*.gz"|wc -l)
count_samp=$(ls dataset|grep -v gz|wc -l)
if [[ "$count_files" -eq 0 ]]; then
    echo -e "\nThere are no sequencing files (*.gz) in the dataset folder!"
    echo -e "Please download or move your paired end files into sample specific subfolders within the dataset folder.\n"
    while true; do
        read -p "Do you wish to download a 3 sample dataset using the downloadToy Snakefile rule? ~1.8 GB of storage required (y/n)" yn
        case $yn in
            [Yy]* ) echo "Running the downloadToy snakefile rule ... " && snakemake downloadToy -j1; break;;
            [Nn]* ) echo "Skipping toy dataset download ... "; break;;
            * ) echo "Please answer yes or no.";;
        esac
    done
elif [[ "$count_samp" -eq 0 && "$count_files" -ne 0 ]]; then 
    echo -e "\nDetected $count_files unorganized files (*.gz) in dataset folder, running organizeData rule ... "
    while true; do
        read -p "Do you wish to organize your samples using the organizeData Snakefile rule? (y/n)" yn
        case $yn in
            [Yy]* ) echo "Running the organizeData snakefile rule ... " && snakemake organizeData -j1; break;;
            [Nn]* ) echo "Skipping toy dataset download ... "; break;;
            * ) echo "Please answer yes or no.";;
        esac
    done
elif [[ "$count_samp" -ne 0 && "$count_files" -ne 0 ]]; then
    echo -e "\nFiles appear to be organized into sample specific subdirectories within the dataset folder."
    echo -e "\nPrinting sample IDs for user verification: "
    ls dataset|grep -v gz
    echo ""
fi

# scratch dir
echo -e "\nPlease remember to set the scratch/ path in the config.yaml file"
echo 'Ideally this path should be set to a job-specific variable that points to a location on your cluster for high I/O operations (e.g. $SCRATCH or $TMPDIR)'
echo "However, it can also be a static directory and metaGEM will create job specific subdirectories automatically."

}

# Run stats task
run_stats() {

echo -e "Checking status of current metaGEM analysis ... \n"

#dataset: count subfolders to determine total number of samples
nsamp=$(ls -d dataset/*/|wc -l)
echo "Raw data: $nsamp samples were identified in the dataset folder ... "

#qfilter: count .json report files
nqfilt=$(find qfiltered -name "*.json"|wc -l)
echo "Quality filtering: $nqfilt / $nsamp samples processed ... "
    
#assembly: count .gz fasta files
nassm=$(find assemblies -name "*.gz"|wc -l)
echo "Assembly: $nassm / $nsamp samples processed ... "
    
#concoct: count *concoct-bins subfolders
nconc=$(find concoct -name "*.concoct-bins"|wc -l)
echo "Binning (CONCOCT): $nconc / $nsamp samples processed ... "
    
#maxbin2: count *maxbin-bins subfolders
nmaxb=$(find maxbin -name "*.maxbin-bins"|wc -l)
echo "Binning (MaxBin2): $nmaxb / $nsamp samples processed ... "
    
#metabat2: count *metabat-bins subfolders
nmetab=$(find metabat -name "*.metabat-bins"|wc -l)
echo "Binning (MetaBAT2): $nmetab / $nsamp samples processed ... "

#metawrap_refine: count subfolders
nmwref=$(ls -d refined_bins/*|wc -l)
echo "Bin refinement: $nmwref / $nsamp samples processed ... "
    
#metawrap_reassemble: count subfolders, also determine total number of final MAGs across samples
nmwrea=$(ls -d reassembled_bins/*|wc -l)
echo "Bin reassembly: $nmwrea / $nsamp samples processed ... "

#taxonomy: count subfolders
ntax=$(ls -d GTDBtk/*|wc -l)
echo "Taxonomy: $ntax / $nsamp samples processed ... "
    
#abundances: count subfolders
nabund=$(ls -d abundance/*|wc -l)
echo "Abundance: $nabund / $nsamp samples processed ... "
    
#models: count subfolders for sample progress and count .xml GEM files for total models generated
ngems=$(find GEMs -name "*xml"|wc -l)
ngemsamp=$(ls -d GEMs/*|wc -l)
echo "GEMs: $ngems models generated from $ngemsamp samples ... "
    
#model reports: count subfolders
nmemo=$(find memote -name "*.gz"|wc -l)
echo "GEM Reports: $nmemo / $ngems models samples ... "
    
#simulations: count .tsv files
nsmet=$(find memote -name "*.gz"|wc -l)
echo -e "GEM Reports: $nsmet / $ngemsamp communities simulated ... \n"

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
    echo -e "\nPlease pay close attention to make sure that your paths are properly configured!"

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

# Submit login node function, note that is only works for rules with no wildcard expansion
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

# Submit local function, similar to submitLogin() but can handle wildcard expanded rules for non-cluster usage
submitLocal() {

    # Parse Snakefile rule all (line 22 of Snakefile) input to match output of desired target rule stored in "$string". Note: Hardcoded line number.
    echo "Parsing Snakefile to target rule: $task ... "
    sed  -i "22s~^.*$~        $string~" Snakefile

    checkParams

    snakeConfig

    echo "Unlocking snakemake ... "
    snakemake --unlock -j 1

    echo -e "\nDry-running snakemake jobs ... "
    snakemake all -n

    while true; do
        read -p "Do you wish to submit this batch of jobs on your local machine? (y/n)" yn
        case $yn in
            [Yy]* ) echo "snakemake all -j 1 -k"|bash; break;;
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
  echo -e "Setting current directory to root in config.yaml file ... \n"
  root=$(pwd)
  sed  -i "2s~/.*$~$root~" config.yaml # hardcoded line for root, change the number 2 if any new lines are added to the start of config.yaml

  # No need to parse snakefile for login node jobs, submit the following locally
  if [ $task == "createFolders" ] || [ $task == "downloadToy" ] || [ $task == "organizeData" ] || [ $task == "qfilterVis" ] || [ $task == "assemblyVis" ] || [ $task == "binningVis" ] || [ $task == "taxonomyVis" ] || [ $task == "abundanceVis" ] || [ $task == "extractProteinBins" ] || [ $task == "extractDnaBins" ] || [ $task == "organizeGEMs" ] || [ $task == "modelVis" ] || [ $task == "interactionVis" ] || [ $task == "growthVis" ] || [ $task == "prepareRoary" ]; then
    submitLogin

  elif [ $task == "check" ]; then
    run_check

  elif [ $task == "stats" ]; then
    run_stats

 # Parse snakefile for cluster/local jobs
  elif [ $task == "fastp" ]; then
    string='expand(config["path"]["root"]+"/"+config["folder"]["qfiltered"]+"/{IDs}/{IDs}_R1.fastq.gz", IDs = IDs)'
    if [ $local == "true" ]; then
        submitLocal
    else
        submitCluster
    fi

  elif [ $task == "megahit" ]; then
    string='expand(config["path"]["root"]+"/"+config["folder"]["assemblies"]+"/{IDs}/contigs.fasta.gz", IDs = IDs)'
    if [ $local == "true" ]; then
        submitLocal
    else
        submitCluster
    fi

  elif [ $task == "crossMapSeries" ]; then
    string='expand(config["path"]["root"]+"/"+config["folder"]["concoct"]+"/{IDs}/cov", IDs = IDs)'
    if [ $local == "true" ]; then
        submitLocal
    else
        submitCluster
    fi

  elif [ $task == "kallistoIndex" ]; then
    string='expand(config["path"]["root"]+"/"+config["folder"]["kallistoIndex"]+"/{focal}/index.kaix", focal = focal)'
    if [ $local == "true" ]; then
        submitLocal
    else
        submitCluster
    fi

  elif [ $task == "crossMapParallel" ]; then
    string='expand(config["path"]["root"]+"/"+config["folder"]["kallisto"]+"/{focal}/{IDs}", focal = focal , IDs = IDs)'
    if [ $local == "true" ]; then
        submitLocal
    else
        submitCluster
    fi

  elif [ $task == "kallisto2concoct" ]; then
    string='expand(config["path"]["root"]+"/"+config["folder"]["concoct"]+"/{focal}/cov/coverage_table.tsv", focal = focal)'
    if [ $local == "true" ]; then
        submitLocal
    else
        submitCluster
    fi

  elif [ $task == "concoct" ]; then
    string='expand(config["path"]["root"]+"/"+config["folder"]["concoct"]+"/{IDs}/{IDs}.concoct-bins", IDs = IDs)'
    if [ $local == "true" ]; then
        submitLocal
    else
        submitCluster
    fi

  elif [ $task == "metabat" ]; then
    string='expand(config["path"]["root"]+"/"+config["folder"]["metabat"]+"/{IDs}/{IDs}.metabat-bins", IDs = IDs)'
    if [ $local == "true" ]; then
        submitLocal
    else
        submitCluster
    fi

  elif [ $task == "maxbin" ]; then
    string='expand(config["path"]["root"]+"/"+config["folder"]["maxbin"]+"/{IDs}/{IDs}.maxbin-bins", IDs = IDs)'
    if [ $local == "true" ]; then
        submitLocal
    else
        submitCluster
    fi

  elif [ $task == "binRefine" ]; then
    string='expand(config["path"]["root"]+"/"+config["folder"]["refined"]+"/{IDs}", IDs = IDs)'
    if [ $local == "true" ]; then
        submitLocal
    else
        submitCluster
    fi

  elif [ $task == "binReassemble" ]; then
    string='expand(config["path"]["root"]+"/"+config["folder"]["reassembled"]+"/{IDs}", IDs = IDs)'
    if [ $local == "true" ]; then
        submitLocal
    else
        submitCluster
    fi

  elif [ $task == "gtdbtk" ]; then
    string='expand(config["path"]["root"]+"/"+config["folder"]["classification"]+"/{IDs}", IDs = IDs)'
    if [ $local == "true" ]; then
        submitLocal
    else
        submitCluster
    fi

  elif [ $task == "abundance" ]; then
    string='expand(config["path"]["root"]+"/"+config["folder"]["abundance"]+"/{IDs}", IDs = IDs)'
    if [ $local == "true" ]; then
        submitLocal
    else
        submitCluster
    fi

  elif [ $task == "carveme" ]; then
    string='expand(config["path"]["root"]+"/"+config["folder"]["GEMs"]+"/{binIDs}.xml", binIDs = binIDs)'
    if [ $local == "true" ]; then
        submitLocal
    else
        submitCluster
    fi

  elif [ $task == "smetana" ]; then
    string='expand(config["path"]["root"]+"/"+config["folder"]["SMETANA"]+"/{IDs}_detailed.tsv", IDs = IDs)'
    if [ $local == "true" ]; then
        submitLocal
    else
        submitCluster
    fi

  elif [ $task == "memote" ]; then
    string='expand(config["path"]["root"]+"/"+config["folder"]["memote"]+"/{gemIDs}", gemIDs = gemIDs)'
    if [ $local == "true" ]; then
        submitLocal
    else
        submitCluster
    fi

  elif [ $task == "grid" ]; then
    string='expand(config["path"]["root"]+"/"+config["folder"]["GRiD"]+"/{IDs}", IDs = IDs)'
    if [ $local == "true" ]; then
        submitLocal
    else
        submitCluster
    fi

  elif [ $task == "prokka" ]; then
    string='expand(config["path"]["root"]+"/"+config["folder"]["pangenome"]+"/prokka/unorganized/{binIDs}", binIDs = binIDs)'
    if [ $local == "true" ]; then
        submitLocal
    else
        submitCluster
    fi

  elif [ $task == "roary" ]; then
    string='expand(config["path"]["root"]+"/"+config["folder"]["pangenome"]+"/roary/{speciesIDs}/", speciesIDs = speciesIDs)'
    if [ $local == "true" ]; then
        submitLocal
    else
        submitCluster
    fi

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
    local=false;
    # Read in options
    while [[ $1 = -?* ]]; do
      case $1 in
        -t|--task) shift; task=${1} ;;
        -j|--nJobs) shift; njobs=${1} ;;
        -c|--nCores) shift; ncores=${1} ;;
        -m|--mem) shift; mem=${1} ;;
        -h|--hours) shift; hours=${1} ;;
        -l|--local) shift; local=true;;
        --endopts) shift; break ;;
        * ) echo "Unknown option(s) provided, please read helpfile ... " && usage && exit 1;;
      esac
      shift
    done
    parse

fi
