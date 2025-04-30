#!/usr/bin/bash
#SBATCH --job-name="02-process_frags"
#SBATCH --time=01:00:00
#SBATCH --output=../../logs/03-chrombpnet/00/%x-%j.out
#SBATCH --partition=akundaje
#SBATCH --cpus-per-task=8
#SBATCH --mem=12G

# load conda environment
eval "$(conda shell.bash hook)"
conda activate chrombpnet

# source configuration variables
source ../config.sh

input_parallel=8
datasets1=$(ls ${sample_frags_dir}/fragments/*.tsv | xargs -n 1 -I {} basename {} .tsv)	
# echo $datasets

datasets2=$( for dataset in ${datasets1[@]}; do

    # check if the out file exists
    # if not, keep the dataset & get list of all unique ones
    done_file="${sample_frags_dir}/pseudorepT/${dataset}.tsv"
    if [[ ! -f "$done_file" ]]; then
        echo "${dataset}.tsv"
    fi

	done | uniq )


echo "@ processing fragments file: ${datasets2}"

parallel -j ${input_parallel} python 02-process_frags.py {} ${sample_frags_dir} ::: ${datasets2}

echo "@ done"
