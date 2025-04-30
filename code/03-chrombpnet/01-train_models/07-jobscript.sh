#!/bin/bash
#SBATCH --output=../../logs/03-chrombpnet/01/07/%x-%j.out
#SBATCH -p akundaje,owners,wjg,sfgf,biochem
#SBATCH -t 01:00:00
#SBATCH -c 4
#SBATCH --mem=20G
#SBATCH -C NO_GPU

celltype=${1}
counts_shaps=${2}
counts_out=${3}
profile_shaps=${4}
profile_out=${5}
peaks_file=${6}
chrom_sizes=${7}

echo -e "@ Using counts shaps: ${counts_shaps} --> ${counts_out}.bw"
echo -e "@ Using profile shaps: ${profile_shaps} --> ${profile_out}.bw"

# source configuration variables
source ../config.sh

# load conda environment
eval "$(conda shell.bash hook)"
conda activate chrombpnet

# the path to the helper script in the ChromBPNet repo
script_loc="${chrombpnet_code}/chrombpnet/evaluation/make_bigwigs/importance_hdf5_to_bigwig.py"

# generate bigwigs if they don't exist yet
if [[ ! -f "${counts_out}.bw" ]]; then
  echo "[$(date +"%m/%d/%Y (%r)")] starting ${celltype} - generating bigwig for counts"
  echo ${counts_shaps}
  echo ${peaks_file}
  echo ${counts_out}
  python3.8 ${script_loc} --hdf5 ${counts_shaps} \
                          --regions ${peaks_file} \
                          --chrom-sizes ${chrom_sizes} \
                          --output-prefix ${counts_out}
fi

echo "[$(date +"%m/%d/%Y (%r)")] done!"

