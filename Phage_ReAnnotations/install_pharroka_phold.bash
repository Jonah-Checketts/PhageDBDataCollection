#!/bin/bash
#Do this first (dont have to rerun it after it has been ran once, provided you are in the conda env pharokka_env)
set -e

# Source conda to allow environment activation in this script
source "$(conda info --base)/etc/profile.d/conda.sh"

PYTHON_VERSION="3.10"
PHAROKKA_VERSION="1.7.5"
PHOLD_VERSION="0.2.0"

echo "python version ${PYTHON_VERSION}"

# new env
if [ ! -f CONDA_READY ]; then
  echo "Creating conda environment 'pharokka_env'"
  conda create -n pharokka_env python=${PYTHON_VERSION} -y
  conda activate pharokka_env
  touch CONDA_READY 
else
  conda activate pharokka_env
fi

# installing pharokka and phold
if [ ! -f PHAROKKA_PHOLD_READY ]; then
  echo "Installing pharokka and phold"
  conda install -y -c conda-forge -c bioconda pip pharokka==${PHAROKKA_VERSION} phold==${PHOLD_VERSION} pytorch
  touch PHAROKKA_PHOLD_READY
fi

echo "Installation complete!"
