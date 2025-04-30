#!/bin/bash
#SBATCH --output=../../logs/03-chrombpnet/01/06/%x-%j.out
#SBATCH -p akundaje,wjg,biochem,sfgf
#SBATCH -t 2-0
#SBATCH --mem=50G
#SBATCH -C NO_GPU
#SBATCH --cpus-per-task=12


peak_shaps=${1}
max_seqlets=${2}
num_leiden=${3}
modisco_output=${4}
report_outdir=${5}
img_suffix_dir=${6}
meme_db=${7}
num_matches=${8}
output_memedb=${9}

# load conda environment
eval "$(conda shell.bash hook)"
conda activate modiscolite

# module load cuda/11.2
# module load cudnn/8.1
module load system
module load libxml2
module load libxslt
module load perl
module load zlib
module load ghostscript
module load cairo

export PATH=$HOME/meme/bin:$HOME/meme/libexec/meme-5.5.5:$PATH

echo "[$(date +"%m/%d/%Y (%r)")] running modisco motifs..."
modisco motifs -i ${peak_shaps} -n ${max_seqlets} -l ${num_leiden} -o ${modisco_output}

echo "[$(date +"%m/%d/%Y (%r)")] running modisco report..."
modisco report -i ${modisco_output} -o ${report_outdir} -s ${img_suffix_dir} -m ${meme_db} -n ${num_matches}

echo "[$(date +"%m/%d/%Y (%r)")] exporting modisco meme..."
modisco meme -i ${modisco_output} -t PFM -o ${output_memedb}


