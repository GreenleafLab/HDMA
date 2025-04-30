#!/bin/bash
#SBATCH --output=../../logs/03-chrombpnet/02/02/%x-%j.out
#SBATCH -p akundaje,wjg,biochem,sfgf
#SBATCH -t 03:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=60G

set -euo pipefail

input=${1}
out_dir=${2}
t=${3}

# source configuration variables
source ../config.sh

# load conda env	
eval "$(conda shell.bash hook)"
conda activate gimme

echo "@ running: gimme cluster ${input} ${out_dir} ${t}"

gimme cluster ${input} ${out_dir} -t ${t} -N 16
