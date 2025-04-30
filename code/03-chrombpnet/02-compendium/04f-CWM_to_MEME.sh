#!/usr/bin/bash

# Purpose: extract CWMs from merged MoDISco patterns in MEME format for downstream
# use in R and for visualization.

#SBATCH --output=../../logs/03-chrombpnet/02/04/04f-cwm_to_meme-%j.out
#SBATCH -p akundaje,wjg,biochem,sfgf
#SBATCH -t 04:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=10G

# fail explicitly for any errors, nonexistent variables, etc
# https://gist.github.com/mohanpedala/1e2ff5661761d3abd0385e8223e16425#set--u
set -euo pipefail

source ../config.sh

# load conda env
eval "$(conda shell.bash hook)"
conda activate modiscolite

module load system
module load cairo

compiled_h5="${modisco_comp_dir}/modisco_compiled.h5"
echo $compiled_h5

python -u 04f-CWM_to_MEME.py --filename $compiled_h5 \
  --datatype CWM \
  --output "${modisco_comp_dir}/modisco_compiled.memedb.txt"

echo "done."
