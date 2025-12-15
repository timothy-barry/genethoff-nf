#!/bin/bash
#SBATCH --job-name=run_nf_pipeline
#SBATCH --output=logs/bioinf_pipeline_%j.log
#SBATCH --error=logs/bioinf_pipeline_%j.log
#SBATCH --time=24:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=4
#SBATCH --partition=short
#SBATCH -c 4

launch_dir="/Users/timbarry/research_code/genethoff-nf/demo/"

nextflow pull timothy-barry/genethoff-nf
nextflow run timothy-barry/genethoff-nf -r nature-revision \
-c nextflow.config \
--conda_env $launch_dir"crisprde-venv.yaml"