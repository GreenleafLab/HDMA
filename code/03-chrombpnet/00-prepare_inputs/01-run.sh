#!/usr/bin/bash
#SBATCH --job-name="01-frags"
#SBATCH --time=02:00:00
#SBATCH --output=../../logs/03-chrombpnet/00/%x-%j.out
#SBATCH --partition=akundaje
#SBATCH --cpus-per-task=4
#SBATCH --mem=110G

# set lib

renv_dir=""
IFS= read -r renv_dir < "../../SJ_RENV_DIR.txt"
export R_LIBS_USER=$renv_dir

# source configuration variables
source ../config.sh

Rscript 01-get_fragments.R $base_dir $sample_frags_dir