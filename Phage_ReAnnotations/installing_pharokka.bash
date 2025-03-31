#!/bin/bash
# Installation script for Pharokka and its database

set -e

# conda environment activation for script
source "$(conda info --base)/etc/profile.d/conda.sh"

PYTHON_VERSION="3.10"
PHAROKKA_VERSION="1.7.5"

echo "Python version ${PYTHON_VERSION}"

# Step 1: Create conda environment if not already created. Need to have it active when doing all annotation stuff.
if [ ! -f CONDA_READY ]; then
  echo "Creating conda environment 'pharokka_env'"
  conda create -n pharokka_env python=${PYTHON_VERSION} -y
  conda activate pharokka_env
  touch CONDA_READY
else
  conda activate pharokka_env
fi

# Step 2: Install pharokka if not already installed
if [ ! -f PHAROKKA_READY ]; then
  echo "Installing pharokka"
  conda install -y -c conda-forge -c bioconda pip pharokka==${PHAROKKA_VERSION} pytorch
  touch PHAROKKA_READY
fi

# Step 3: Download pharokka database (only needs to be done once, is a large database)
if [ ! -f PHAROKKA_DB_READY ]; then
  echo "Downloading pharokka database. This will take a few minutes. Please be patient :)"
  time install_databases.py -o pharokka_db
  touch PHAROKKA_DB_READY
fi

echo "Pharokka installation and database setup complete!"
