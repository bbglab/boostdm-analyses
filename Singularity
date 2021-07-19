Bootstrap: docker

From: ubuntu

%help
    Singularity image to run boostdm-analyses notebooks with all dependencies
    GitHub: https://github.com/bbglab/boostdm-analyses

%setup
    mkdir ${SINGULARITY_ROOTFS}/notebook

%post
    apt-get update
    apt-get -y --no-install-recommends install python3 python3-dev python3-pip python3-tqdm python3-numpy python3-pandas 
    /usr/bin/pip3 install ipython jupyterlab
