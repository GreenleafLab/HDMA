#!/usr/bin/bash
#SBATCH --job-name="05-format_peaks"
#SBATCH --time=01:00:00
#SBATCH --output=../../logs/03-chrombpnet/00/%x-%j.out
#SBATCH --partition=akundaje
#SBATCH --cpus-per-task=8
#SBATCH --mem=5G

# load conda environment
eval "$(conda shell.bash hook)"
conda activate chrombpnet

# source configuration variables
source ../config.sh

out_dir=$chrombpnet_peaks_dir
input_parallel=8
datasets=$(ls $peaks_dir/*__peaks_overlap_filtered.narrowPeak | xargs -n 1 -I {} basename {} __peaks_overlap_filtered.narrowPeak)	

echo "@ formatting peaks for clusters: ${datasets}"

# DEBUG:
# python 05-format_peaks.py Adrenal_c0 ${peaks_dir} ${out_dir}

parallel -j ${input_parallel} python 05-format_peaks.py {} ${peaks_dir} ${out_dir} ::: ${datasets}

echo "@ done"
