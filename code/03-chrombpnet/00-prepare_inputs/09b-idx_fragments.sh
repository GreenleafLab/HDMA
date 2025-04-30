#!/usr/bin/bash
#SBATCH --job-name="10b-idx_frags"
#SBATCH --time=06:00:00
#SBATCH --output=../../logs/03-chrombpnet/00/%x-%j.out
#SBATCH --partition=akundaje,wjg,sfgf,biochem
#SBATCH --cpus-per-task=8
#SBATCH --mem=10G

# Compress and sort peaks for space saving.

# load conda environment
eval "$(conda shell.bash hook)"
conda activate chrombpnet

# source configuration variables
source ../config.sh

module load biology samtools

input_parallel=8

fragsets=$(ls ${cluster_frags_dir}/fragments/*.tsv )
echo ${fragsets[@]}

idx_frag () {

    fragset=$1
    bgzip -c ${fragset} > ${fragset}.gz
    tabix -p bed ${fragset}.gz

}
export -f idx_frag

parallel -j ${input_parallel} idx_frag {} ::: ${fragsets}
