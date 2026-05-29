#!/bin/bash
#SBATCH --job-name=classify-oil-red-o-image
#SBATCH -p kibr,owners
#SBATCH --output=logs/job_%A_%a.out
#SBATCH --error=logs/job_%A_%a.err
#SBATCH --array=0-7 # change to the number of items in the array (may need updating)
#SBATCH --time=48:00:00
#SBATCH --mem=200000 
#SBATCH -c 5
#SBATCH -n 1
#SBATCH -N 1

ml python/3.12.1

source .venv/bin/activate
cd /oak/stanford/projects/kibr/Reorganizing/Projects/James/lipid-droplet-pipeline
python3 /oak/stanford/projects/kibr/Reorganizing/Projects/James/lipid-droplet-pipeline/sbatch/random-forest-image-classifier.py $SLURM_ARRAY_TASK_ID