#!/usr/bin/bash
#SBATCH --job-name="03b-sync"
#SBATCH --time=04:00:00
#SBATCH --output=../../logs/03-chrombpnet/00/%x-%j.out
#SBATCH --partition=akundaje
#SBATCH --cpus-per-task=1
#SBATCH --mem=30G

# source configuration variables
source ../config.sh

rsync -rv ${cluster_frags_dir%/} \
	"${base_dir}/00-inputs"
