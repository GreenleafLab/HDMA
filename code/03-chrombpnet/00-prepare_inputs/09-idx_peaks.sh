#!/usr/bin/bash
#SBATCH --job-name="09-idx_peaks"
#SBATCH --time=02:00:00
#SBATCH --output=../../logs/03-chrombpnet/00/%x-%j.out
#SBATCH --partition=akundaje,wjg,sfgf,biochem
#SBATCH --cpus-per-task=8
#SBATCH --mem=10G

# Compress and sort peaks for viewing in the WashU browser.

# load conda environment
eval "$(conda shell.bash hook)"
conda activate chrombpnet

# source configuration variables
source ../config.sh

module load biology samtools

input_parallel=8

peaksets=$(ls ${chrombpnet_peaks_dir}/*.narrowPeak )
echo ${peaksets[@]}

idx_pk () {

    peakset=$1
    bgzip -c ${peakset} > ${peakset}.gz
    tabix -p bed ${peakset}.gz

}
export -f idx_pk

parallel -j ${input_parallel} idx_pk {} ::: ${peaksets}
