#!/usr/bin/bash
#SBATCH --job-name="03d-sync"
#SBATCH --time=06:00:00
#SBATCH --output=../../logs/03-chrombpnet/00/%x-%j.out
#SBATCH --partition=akundaje
#SBATCH --cpus-per-task=1
#SBATCH --mem=60G

# source configuration variables
source ../config.sh

rsync -rv --ignore-existing \
	"${base_dir}/00-inputs/cluster_fragments/" \
	$cluster_frags_dir
