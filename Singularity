Bootstrap: docker
From: continuumio/miniconda3

%labels
    DESCRIPTION Singularity image containing the conda environment for metaBAGpipes

%environment
    PATH=/opt/conda/envs/metabagpipes/bin:$PATH
    export PATH

%files
    metaBAGpipes_env.yml /

%post
    /opt/conda/bin/conda env create -f /metaBAGpipes_env.yml
    /opt/conda/bin/conda clean -a
    mkdir -p /c3se
    mkdir -p /local
    mkdir -p /apps
    mkdir -p /usr/share/lmod/lmod
    mkdir -p /var/hasplm
    mkdir -p /var/opt/thinlinc
    mkdir -p /usr/lib64
    touch /usr/lib64/libdlfaker.so
    touch /usr/lib64/libvglfaker.so
    touch /usr/bin/nvidia-smi

