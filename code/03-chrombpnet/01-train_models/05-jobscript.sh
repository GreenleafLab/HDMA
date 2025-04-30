#!/bin/bash
#SBATCH --output=../../logs/03-chrombpnet/01/05/%x-%j.out
#SBATCH -p akundaje,sfgf,wjg,biochem
#SBATCH -t 06:00:00
#SBATCH -c 1
#SBATCH --mem=60G
#SBATCH -C NO_GPU

# load conda environment
eval "$(conda shell.bash hook)"
conda activate chrombpnet

# get vars
dir=${1}
folds=${2}

echo ${dir}
echo ${folds}

# run for one cell type and one type of Shap scores
if [[ ! -f "${dir}/average_shaps.counts.h5" ]]; then
  echo "@ averaging COUNTS..."
  python 05-average_deepshaps.py ${dir} ${folds} counts
fi
