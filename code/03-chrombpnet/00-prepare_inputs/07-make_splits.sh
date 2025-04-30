#!/bin/bash
#SBATCH --job-name="make_splits"
#SBATCH --time=01:00:00
#SBATCH --output=../../logs/03-chrombpnet/00/%x-%j.out
#SBATCH --partition=akundaje
#SBATCH --cpus-per-task=1
#SBATCH --mem=5G
#SBATCH -C NO_GPU

# This file creates the ChromBPNet chromosomal folds files specifying train/test/validation
# folds, using the splits in https://zenodo.org/records/7445373

# load conda environment
eval "$(conda shell.bash hook)"
conda activate chrombpnet

# get configuration variables
source ../config.sh

chrombpnet prep splits -c "${refs}/GRCh38_EBV_sorted_standard.chrom.sizes.tsv" --test-chroms chr1 chr3 chr6 --valid-chroms chr8 chr20 -op "${split_dir}/fold_0"
chrombpnet prep splits -c "${refs}/GRCh38_EBV_sorted_standard.chrom.sizes.tsv" --test-chroms chr2 chr8 chr9 chr16 --valid-chroms chr12 chr17 -op "${split_dir}/fold_1"
chrombpnet prep splits -c "${refs}/GRCh38_EBV_sorted_standard.chrom.sizes.tsv" --test-chroms chr4 chr11 chr12 chr15 chrY --valid-chroms chr22 chr7 -op "${split_dir}/fold_2"
chrombpnet prep splits -c "${refs}/GRCh38_EBV_sorted_standard.chrom.sizes.tsv" --test-chroms chr5 chr10 chr14 chr18 chr20 chr22 --valid-chroms chr6 chr21 -op "${split_dir}/fold_3"
chrombpnet prep splits -c "${refs}/GRCh38_EBV_sorted_standard.chrom.sizes.tsv" --test-chroms chr7 chr13 chr17 chr19 chr21 chrX --valid-chroms chr10 chr18 -op "${split_dir}/fold_4"
