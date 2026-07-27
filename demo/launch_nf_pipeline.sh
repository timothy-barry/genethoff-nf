#!/bin/bash
#SBATCH --job-name=run_nf_pipeline
#SBATCH --output=logs/bioinf_pipeline_%j.log
#SBATCH --error=logs/bioinf_pipeline_%j.log
#SBATCH --time=24:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=4
#SBATCH --partition=short
#SBATCH -c 4

nextflow run "/Users/timbarry/research_code/genethoff-nf/main.nf" -c nextflow.config --pool_samples true