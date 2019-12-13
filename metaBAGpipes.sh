#!/bin/bash

# Helpfile function

usage() {
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
                            moveBins
                            carveme
                            organizeGems
                            smetana
                            memote   
                            grid

                        VISUALIZATION (in development)
                            qcVis
                            assemblyVis
                            binningVis
                            taxonomyVis
                            modelVis
                            interactionVis
                            growthVis

  -j, --nJobs       Specify number of jobs to run in parallel
  -c, --nCores      Specify number of cores per job
  -m, --mem         Specify memory in GB required for job
  -h, --help        Display this help and exit

Example: bash metaBAGpipes.sh -t createFolders -j 1 -c 1

"
}

# Parse and submit function

function parse() {

  # Set root folder

  root=$(pwd)
  sed  -i "2s~/.*$~$root~" config.yaml
  echo "Current directory set to root."

  echo "Parsing Snakefile and cluster_config.json to target rule: $task."

  if [ $task == "createFolders" ]; then
    echo "No need to parse Snakefile for rule $task."
    echo "Running snakemake on login node."
    snakemake --unlock
    snakemake $task -n
    while true; do
        read -p "Do you wish to submit this job? (y/n)" yn
        case $yn in
            [Yy]* ) snakemake $task ; break;;
            [Nn]* ) exit;;
            * ) echo "Please answer yes or no.";;
        esac
    done

  elif [ $task == "downloadToy" ]; then
    echo "No need to parse Snakefile for rule $task."
    echo "Running snakemake on login node."
    snakemake --unlock
    snakemake $task -n
    while true; do
        read -p "Do you wish to submit this job? (y/n)" yn
        case $yn in
            [Yy]* ) snakemake $task ; break;;
            [Nn]* ) exit;;
            * ) echo "Please answer yes or no.";;
        esac
    done

  elif [ $task == "organizeData" ]; then
    echo "No need to parse Snakefile for rule $task."
    echo "Running snakemake on login node."
    snakemake --unlock
    snakemake $task -n
    while true; do
        read -p "Do you wish to submit this job? (y/n)" yn
        case $yn in
            [Yy]* ) snakemake $task ; break;;
            [Nn]* ) exit;;
            * ) echo "Please answer yes or no.";;
        esac
    done

  elif [ $task == "fastp" ]; then
    string='expand(config["path"]["root"]+"/"+config["folder"]["qfiltered"]+"/{IDs}/{IDs}_1.fastq.gz", IDs = IDs)'
    sed  -i "15s~^.*$~        $string~" Snakefile
    sed -i "5s/:.*$/: $ncores,/" cluster_config.json
    echo "Submitting $njobs jobs with $ncores cores and $mem memory each."
    echo "Note: Errors in snakemake submission of jobs may likely be due to wildcard missassignment caused by any non-sample ID folder/file in the dataset folder."
    snakemake --unlock
    snakemake all -j $njobs -n -k --cluster-config cluster_config.json -c "sbatch -A {cluster.account} -t {cluster.time} -n {cluster.n} --ntasks {cluster.tasks} --cpus-per-task {cluster.n} --output {cluster.output}"
    while true; do
        read -p "Do you wish to submit this batch of jobs? (y/n)" yn
        case $yn in
            [Yy]* ) echo "nohup snakemake all -j $njobs -k --cluster-config cluster_config.json -c 'sbatch -A {cluster.account} -t {cluster.time} -n {cluster.n} --ntasks {cluster.tasks} --cpus-per-task {cluster.n} --output {cluster.output}' &"|bash; break;;
            [Nn]* ) exit;;
            * ) echo "Please answer yes or no.";;
        esac
    done

  elif [ $task == "megahit" ]; then
    string='expand(config["path"]["root"]+"/"+config["folder"]["assemblies"]+"/{IDs}/contigs.fasta.gz", IDs = IDs)'
    sed  -i "15s~^.*$~        $string~" Snakefile
    sed -i "5s/:.*$/: $ncores,/" cluster_config.json
    sed -i "7s/:.*$/: $(echo $mem)G,/" cluster_config.json
    echo "Submitting $njobs jobs with $ncores cores and $mem memory each."
    echo "Note: Errors in snakemake submission of jobs may likely be due to wildcard missassignment caused by any non-sample ID folder/file in the dataset folder."
    snakemake --unlock
    snakemake all -j $njobs -n -k --cluster-config cluster_config.json -c "sbatch -A {cluster.account} -t {cluster.time} --mem {cluster.mem} -n {cluster.n} --ntasks {cluster.tasks} --cpus-per-task {cluster.n} --output {cluster.output}"
    while true; do
        read -p "Do you wish to submit this batch of jobs? (y/n)" yn
        case $yn in
            [Yy]* ) echo "nohup snakemake all -j $njobs -k --cluster-config cluster_config.json -c 'sbatch -A {cluster.account} -t {cluster.time} --mem {cluster.mem} -n {cluster.n} --ntasks {cluster.tasks} --cpus-per-task {cluster.n} --output {cluster.output}' &"|bash; break;;
            [Nn]* ) exit;;
            * ) echo "Please answer yes or no.";;
        esac
    done

    elif [ $task == "assemblyVis" ]; then
    echo "No need to parse Snakefile for rule $task."
    echo "Running snakemake on login node."
    snakemake --unlock
    snakemake $task -n
    while true; do
        read -p "Do you wish to submit this job? (y/n)" yn
        case $yn in
            [Yy]* ) snakemake $task ; break;;
            [Nn]* ) exit;;
            * ) echo "Please answer yes or no.";;
        esac
    done

  elif [ $task == "kallisto" ]; then
    string='expand(config["path"]["root"]+"/"+config["folder"]["concoctInput"]+"/{IDs}_concoct_inputtableR.tsv", IDs = IDs)'
    sed  -i "15s~^.*$~        $string~" Snakefile
    sed -i "5s/:.*$/: $ncores,/" cluster_config.json
    sed -i "7s/:.*$/: $(echo $mem)G,/" cluster_config.json
    echo "Submitting $njobs jobs with $ncores cores and $mem memory each."
    snakemake --unlock
    snakemake all -j $njobs -n -k --cluster-config cluster_config.json -c "sbatch -A {cluster.account} -t {cluster.time} -n {cluster.n} --ntasks {cluster.tasks} --cpus-per-task {cluster.n} --output {cluster.output}"
        while true; do
        read -p "Do you wish to submit this batch of jobs? (y/n)" yn
        case $yn in
            [Yy]* ) echo "nohup snakemake all -j $njobs -k --cluster-config cluster_config.json -c 'sbatch -A {cluster.account} -t {cluster.time} --mem {cluster.mem} -n {cluster.n} --ntasks {cluster.tasks} --cpus-per-task {cluster.n} --output {cluster.output}' &"|bash; break;;
            [Nn]* ) exit;;
            * ) echo "Please answer yes or no.";;
        esac
    done

  elif [ $task == "concoct" ]; then
    string='expand(config["path"]["root"]+"/"+config["folder"]["concoctOutput"]+"/{IDs}/{IDs}.concoct-bins", IDs = IDs)'
    sed  -i "15s~^.*$~        $string~" Snakefile
    sed -i "5s/:.*$/: $ncores,/" cluster_config.json
    echo "Submitting $njobs jobs with $ncores cores and $mem memory each."
    snakemake --unlock
    snakemake all -j $njobs -n -k --cluster-config cluster_config.json -c "sbatch -A {cluster.account} -t {cluster.time} -n {cluster.n} --ntasks {cluster.tasks} --cpus-per-task {cluster.n} --output {cluster.output}"
        while true; do
        read -p "Do you wish to submit this batch of jobs? (y/n)" yn
        case $yn in
            [Yy]* ) echo "nohup snakemake all -j $njobs -k --cluster-config cluster_config.json -c 'sbatch -A {cluster.account} -t {cluster.time} -n {cluster.n} --ntasks {cluster.tasks} --cpus-per-task {cluster.n} --output {cluster.output}' &"|bash; break;;
            [Nn]* ) exit;;
            * ) echo "Please answer yes or no.";;
        esac
    done

  elif [ $task == "metabat" ]; then
    string='expand(config["path"]["root"]+"/"+config["folder"]["metabat"]+"/{IDs}/{IDs}.metabat-bins", IDs = IDs)'
    sed  -i "15s~^.*$~        $string~" Snakefile
    sed -i "5s/:.*$/: $ncores,/" cluster_config.json
    echo "Submitting $njobs jobs with $ncores cores and $mem memory each."
    snakemake --unlock
    snakemake all -j $njobs -n -k --cluster-config cluster_config.json -c "sbatch -A {cluster.account} -t {cluster.time} -n {cluster.n} --ntasks {cluster.tasks} --cpus-per-task {cluster.n} --output {cluster.output}"
        while true; do
        read -p "Do you wish to submit this batch of jobs?" yn
        case $yn in
            [Yy]* ) echo "nohup snakemake all -j $njobs -k --cluster-config cluster_config.json -c 'sbatch -A {cluster.account} -t {cluster.time} -n {cluster.n} --ntasks {cluster.tasks} --cpus-per-task {cluster.n} --output {cluster.output}' &"|bash; break;;
            [Nn]* ) exit;;
            * ) echo "Please answer yes or no.";;
        esac
    done

  elif [ $task == "maxbin" ]; then
    string='expand(config["path"]["root"]+"/"+config["folder"]["maxbin"]+"/{IDs}/{IDs}.maxbin-bins", IDs = IDs)'
    sed  -i "15s~^.*$~        $string~" Snakefile
    sed -i "5s/:.*$/: $ncores,/" cluster_config.json
    echo "Submitting $njobs jobs with $ncores cores and $mem memory each."
    snakemake --unlock
    snakemake all -j $njobs -n -k --cluster-config cluster_config.json -c "sbatch -A {cluster.account} -t {cluster.time} -n {cluster.n} --ntasks {cluster.tasks} --cpus-per-task {cluster.n} --output {cluster.output}"
        while true; do
        read -p "Do you wish to submit this batch of jobs? (y/n)" yn
        case $yn in
            [Yy]* ) echo "nohup snakemake all -j $njobs -k --cluster-config cluster_config.json -c 'sbatch -A {cluster.account} -t {cluster.time} -n {cluster.n} --ntasks {cluster.tasks} --cpus-per-task {cluster.n} --output {cluster.output}' &"|bash; break;;
            [Nn]* ) exit;;
            * ) echo "Please answer yes or no.";;
        esac
    done

  elif [ $task == "binRefine" ]; then
    string='expand(config["path"]["root"]+"/"+config["folder"]["refined"]+"/{IDs}", IDs = IDs)'
    sed  -i "15s~^.*$~        $string~" Snakefile
    sed -i "5s/:.*$/: $ncores,/" cluster_config.json
    echo "Submitting $njobs jobs with $ncores cores and $mem memory each."
    snakemake --unlock
    snakemake all -j $njobs -n -k --cluster-config cluster_config.json -c "sbatch -A {cluster.account} -t {cluster.time} -n {cluster.n} --ntasks {cluster.tasks} --cpus-per-task {cluster.n} --output {cluster.output}"
        while true; do
        read -p "Do you wish to submit this batch of jobs? (y/n)" yn
        case $yn in
            [Yy]* ) echo "nohup snakemake all -j $njobs -k --cluster-config cluster_config.json -c 'sbatch -A {cluster.account} -t {cluster.time} -n {cluster.n} --ntasks {cluster.tasks} --cpus-per-task {cluster.n} --output {cluster.output}' &"|bash; break;;
            [Nn]* ) exit;;
            * ) echo "Please answer yes or no.";;
        esac
    done

  elif [ $task == "binReassemble" ]; then
    string='expand(config["path"]["root"]+"/"+config["folder"]["reassembled"]+"/{IDs}", IDs = IDs)'
    sed  -i "15s~^.*$~        $string~" Snakefile
    sed -i "5s/:.*$/: $ncores,/" cluster_config.json
    echo "Submitting $njobs jobs with $ncores cores each."
    snakemake --unlock
    snakemake all -j $njobs -n -k --cluster-config cluster_config.json -c "sbatch -A {cluster.account} -t {cluster.time} -n {cluster.n} --ntasks {cluster.tasks} --cpus-per-task {cluster.n} --output {cluster.output}"
        while true; do
        read -p "Do you wish to submit this batch of jobs? (y/n)" yn
        case $yn in
            [Yy]* ) echo "nohup snakemake all -j $njobs -k --cluster-config cluster_config.json -c 'sbatch -A {cluster.account} -t {cluster.time} -n {cluster.n} --ntasks {cluster.tasks} --cpus-per-task {cluster.n} --output {cluster.output}' &"|bash; break;;
            [Nn]* ) exit;;
            * ) echo "Please answer yes or no.";;
        esac
    done

    elif [ $task == "binningVis" ]; then
    echo "No need to parse Snakefile for rule $task."
    echo "Running snakemake on login node."
    snakemake --unlock
    snakemake $task -n
    while true; do
        read -p "Do you wish to submit this job? (y/n)" yn
        case $yn in
            [Yy]* ) snakemake $task ; break;;
            [Nn]* ) exit;;
            * ) echo "Please answer yes or no.";;
        esac
    done

  elif [ $task == "classifyGenomes" ]; then
    string='expand(config["path"]["root"]+"/"+config["folder"]["classification"]+"/{IDs}", IDs = IDs)'
    sed  -i "15s~^.*$~        $string~" Snakefile
    sed -i "5s/:.*$/: $ncores,/" cluster_config.json
    echo "Submitting $njobs jobs with $ncores cores and $mem memory each."
    snakemake --unlock
    snakemake all -j $njobs -n -k --cluster-config cluster_config.json -c "sbatch -A {cluster.account} -t {cluster.time} -n {cluster.n} --ntasks {cluster.tasks} --cpus-per-task {cluster.n} --output {cluster.output}"
        while true; do
        read -p "Do you wish to submit this batch of jobs? (y/n)" yn
        case $yn in
            [Yy]* ) echo "nohup snakemake all -j $njobs -k --cluster-config cluster_config.json -c 'sbatch -A {cluster.account} -t {cluster.time} -n {cluster.n} --ntasks {cluster.tasks} --cpus-per-task {cluster.n} --output {cluster.output}' &"|bash; break;;
            [Nn]* ) exit;;
            * ) echo "Please answer yes or no.";;
        esac
    done

  elif [ $task == "abundance" ]; then
    string='expand(config["path"]["root"]+"/"+config["folder"]["abundance"]+"/{IDs}", IDs = IDs)'
    sed  -i "15s~^.*$~        $string~" Snakefile
    sed -i "5s/:.*$/: $ncores,/" cluster_config.json
    echo "Submitting $njobs jobs with $ncores cores and $mem memory each."
    snakemake --unlock
    snakemake all -j $njobs -n -k --cluster-config cluster_config.json -c "sbatch -A {cluster.account} -t {cluster.time} -n {cluster.n} --ntasks {cluster.tasks} --cpus-per-task {cluster.n} --output {cluster.output}"
        while true; do
        read -p "Do you wish to submit this batch of jobs? (y/n)" yn
        case $yn in
            [Yy]* ) echo "nohup snakemake all -j $njobs -k --cluster-config cluster_config.json -c 'sbatch -A {cluster.account} -t {cluster.time} -n {cluster.n} --ntasks {cluster.tasks} --cpus-per-task {cluster.n} --output {cluster.output}' &"|bash; break;;
            [Nn]* ) exit;;
            * ) echo "Please answer yes or no.";;
        esac
    done

  elif [ $task == "moveBins" ]; then
    echo "No need to parse Snakefile for rule $task."
    echo "Running snakemake on login node."
    snakemake --unlock
    snakemake $task -n
    while true; do
        read -p "Do you wish to submit this batch of jobs? (y/n)" yn
        case $yn in
            [Yy]* ) snakemake $task ; break;;
            [Nn]* ) exit;;
            * ) echo "Please answer yes or no.";;
        esac
    done

  elif [ $task == "carveme" ]; then
    string='expand(config["path"]["root"]+"/"+config["folder"]["GEMs"]+"/{binIDs}.xml", binIDs = binIDs)'
    sed  -i "15s~^.*$~        $string~" Snakefile
    sed -i "5s/:.*$/: $ncores,/" cluster_config.json
    echo "Submitting $njobs jobs with $ncores cores and $mem memory each."
    snakemake --unlock
    snakemake all -j $njobs -n -k --cluster-config cluster_config.json -c "sbatch -A {cluster.account} -t {cluster.time} -n {cluster.n} --ntasks {cluster.tasks} --cpus-per-task {cluster.n} --output {cluster.output}"
        while true; do
        read -p "Do you wish to submit this batch of jobs? (y/n)" yn
        case $yn in
            [Yy]* ) echo "nohup snakemake all -j $njobs -k --cluster-config cluster_config.json -c 'sbatch -A {cluster.account} -t {cluster.time} -n {cluster.n} --ntasks {cluster.tasks} --cpus-per-task {cluster.n} --output {cluster.output}' &"|bash; break;;
            [Nn]* ) exit;;
            * ) echo "Please answer yes or no.";;
        esac
    done

  elif [ $task == "organizeGems" ]; then
    echo "No need to parse Snakefile for rule $task."
    echo "Running snakemake on login node."
    snakemake --unlock
    snakemake $task -n
    while true; do
        read -p "Do you wish to submit this job? (y/n)" yn
        case $yn in
            [Yy]* ) snakemake $task ; break;;
            [Nn]* ) exit;;
            * ) echo "Please answer yes or no.";;
        esac
    done

    elif [ $task == "modelVis" ]; then
    echo "No need to parse Snakefile for rule $task."
    echo "Running snakemake on login node."
    snakemake --unlock
    snakemake $task -n
    while true; do
        read -p "Do you wish to submit this job? (y/n)" yn
        case $yn in
            [Yy]* ) snakemake $task ; break;;
            [Nn]* ) exit;;
            * ) echo "Please answer yes or no.";;
        esac
    done

  elif [ $task == "smetana" ]; then
    string='expand(config["path"]["root"]+"/"+config["folder"]["SMETANA"]+"/{IDs}_detailed.tsv", IDs = IDs)'
    sed  -i "15s~^.*$~        $string~" Snakefile
    sed -i "5s/:.*$/: $ncores,/" cluster_config.json
    echo "Submitting $njobs jobs with $ncores cores and $mem memory each."
    snakemake --unlock
    snakemake all -j $njobs -n -k --cluster-config cluster_config.json -c "sbatch -A {cluster.account} -t {cluster.time} -n {cluster.n} --ntasks {cluster.tasks} --cpus-per-task {cluster.n} --output {cluster.output}"
        while true; do
        read -p "Do you wish to submit this batch of jobs? (y/n)" yn
        case $yn in
            [Yy]* ) echo "nohup snakemake all -j $njobs -k --cluster-config cluster_config.json -c 'sbatch -A {cluster.account} -t {cluster.time} -n {cluster.n} --ntasks {cluster.tasks} --cpus-per-task {cluster.n} --output {cluster.output}' &"|bash; break;;
            [Nn]* ) exit;;
            * ) echo "Please answer yes or no.";;
        esac
    done

  elif [ $task == "memote" ]; then
    string='expand(config["path"]["root"]+"/"+config["folder"]["memote"]+"/{IDs}", IDs = IDs)'
    sed  -i "15s~^.*$~        $string~" Snakefile
    sed -i "5s/:.*$/: $ncores,/" cluster_config.json
    echo "Submitting $njobs jobs with $ncores cores and $mem memory each."
    snakemake --unlock
    snakemake all -j $njobs -n -k --cluster-config cluster_config.json -c "sbatch -A {cluster.account} -t {cluster.time} -n {cluster.n} --ntasks {cluster.tasks} --cpus-per-task {cluster.n} --output {cluster.output}"
        while true; do
        read -p "Do you wish to submit this batch of jobs? (y/n)" yn
        case $yn in
            [Yy]* ) echo "nohup snakemake all -j $njobs -k --cluster-config cluster_config.json -c 'sbatch -A {cluster.account} -t {cluster.time} -n {cluster.n} --ntasks {cluster.tasks} --cpus-per-task {cluster.n} --output {cluster.output}' &"|bash; break;;
            [Nn]* ) exit;;
            * ) echo "Please answer yes or no.";;
        esac
    done

  elif [ $task == "grid" ]; then
    string='expand(config["path"]["root"]+"/"+config["folder"]["GRiD"]+"/{IDs}", IDs = IDs)'
    sed  -i "15s~^.*$~        $string~" Snakefile
    sed -i "5s/:.*$/: $ncores,/" cluster_config.json
    echo "Submitting $njobs jobs with $ncores cores and $mem memory each."
    snakemake --unlock
    snakemake all -j $njobs -n -k --cluster-config cluster_config.json -c "sbatch -A {cluster.account} -t {cluster.time} -n {cluster.n} --ntasks {cluster.tasks} --cpus-per-task {cluster.n} --output {cluster.output}"
        while true; do
        read -p "Do you wish to submit this batch of jobs? (y/n)" yn
        case $yn in
            [Yy]* ) echo "nohup snakemake all -j $njobs -k --cluster-config cluster_config.json -c 'sbatch -A {cluster.account} -t {cluster.time} -n {cluster.n} --ntasks {cluster.tasks} --cpus-per-task {cluster.n} --output {cluster.output}' &"|bash; break;;
            [Nn]* ) exit;;
            * ) echo "Please answer yes or no.";;
        esac
    done

  else
    echo "Task not recognized."
    usage
  fi

}

# Check if arguments are provided, if so run parse function

if [ $# -eq 0 ]; then
    echo "No arguments provided ... "
    usage
else
    # Read in options
    while [[ $1 = -?* ]]; do
      case $1 in
        -h|--help) usage >&2;;
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
