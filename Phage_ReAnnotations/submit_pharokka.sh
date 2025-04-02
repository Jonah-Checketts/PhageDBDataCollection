#!/bin/bash
#SBATCH --job-name=pharokka_job       
#SBATCH --output=pharokka_job.out #can view this with tail -f pharokka_job.out
#SBATCH --error=pharokka_job.err
#SBATCH --time=04:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G

# Load/activate whatever environment you need for Pharokka
module load pharokka    # (Or remove if not needed)
mamba activate pharokka_env

# Forward command-line args from the "sbatch" command to run_clusters.bash:
bash run_clusters.bash "$@"
# Clusters that we have existing data for: BE, BU, AU, DJ, CZ, BI, AK, AS, P, BD
#check jobs with squeue -u $USER
#check job history with sacct -u $USER