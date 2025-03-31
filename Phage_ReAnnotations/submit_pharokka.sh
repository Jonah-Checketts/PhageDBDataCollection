#!/bin/bash
#SBATCH --job-name=pharokka_job       
#SBATCH --output=pharokka_job.out      #can view this with tail -f pharokka_job.out
#SBATCH --error=pharokka_job.err
#SBATCH --time=04:00:00                #max time
#SBATCH --ntasks=2                     
#SBATCH --cpus-per-task=8              
#SBATCH --mem=32G         
       

# module load pharokka
#Make sure you are in the correct mamba/conda environment
mamba activate pharokka_env

bash run_pharokka.bash

#Run command: sbatch submit_pharokka.sh
#check jobs with squeue -u $USER
#check job history with sacct -u $USER