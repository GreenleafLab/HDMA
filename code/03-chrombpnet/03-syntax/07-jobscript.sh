#!/usr/bin/bash
#SBATCH --output=../../logs/03-chrombpnet/03/07/%x-%j.out
#SBATCH -p akundaje,wjg,sfgf,biochem
#SBATCH --time=02:00:00
#SBATCH -c 1
#SBATCH --mem=10G

# load conda environment
eval "$(conda shell.bash hook)"
conda activate tangermeme

source ../config.sh

# get all arguments
clusters=("$@")

python -u 07-motif_cooccurrence.py --clusters ${clusters[@]}