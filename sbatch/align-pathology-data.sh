#!/bin/bash
#SBATCH --job-name=align-pathology-data
#SBATCH -p kibr,owners
#SBATCH --output=logs/job_%A_%a.out
#SBATCH --error=logs/job_%A_%a.err
#SBATCH --array=0-11
#SBATCH --time=00:40:00
#SBATCH --mem=20000 
#SBATCH -c 10
#SBATCH -n 1
#SBATCH -N 1

ml python/3.12.1

STRINGS=("04-06" "04-44" "05-27" "10-46" "13-54" "13-69" "14-02" "14-20" "15-27" "18-20" "18-75" "99-15")
CURRENT_STRING=${STRINGS[$SLURM_ARRAY_TASK_ID]}

echo "Processing string: $CURRENT_STRING"

cd /oak/stanford/projects/kibr/Reorganizing/Projects/James/lipid-droplet-pipeline
source .venv/bin/activate
python3 /oak/stanford/projects/kibr/Reorganizing/Projects/James/lipid-droplet-pipeline/sbatch/align-pathology-data.py "$CURRENT_STRING"