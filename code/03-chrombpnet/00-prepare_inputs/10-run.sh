#!/usr/bin/bash
#SBATCH --job-name="10-subset_background"
#SBATCH --time=03:00:00
#SBATCH --output=../../logs/03-chrombpnet/00/%x-%j.out
#SBATCH --partition=akundaje,wjg,owners,sfgf,biochem
#SBATCH --cpus-per-task=8
#SBATCH --mem=5G

# load conda environment
eval "$(conda shell.bash hook)"
conda activate chrombpnet

# source configuration variables
source ../config.sh

input_parallel=8

datasets=$(ls $peaks_dir/*__peaks_overlap_filtered.narrowPeak | xargs -n 1 -I {} basename {} __peaks_overlap_filtered.narrowPeak)	

echo "@ formatting peaks for clusters: ${datasets}"

# DEBUG:
for dataset in ${datasets[@]}; do

  python 13-subset_backgrounds.py --cluster $dataset --negatives-dir ${negatives_dir} --output-dir ${negatives_subset_dir}
  
done

echo "@ done"
