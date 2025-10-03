#!/bin/bash
#SBATCH --job-name=12-annotate_hits
#SBATCH --output=../../logs/03-chrombpnet/00/12/%x-%j.out
#SBATCH --partition=akundaje,wjg,biochem,sfgf
#SBATCH --mem-per-cpu=16G
#SBATCH --cpus-per-task=6
#SBATCH --time=01:00:00

# Purpose: runner script to annotate peaks.
#
# Usage:
# $ organs=( Adrenal Brain Eye Heart Liver Lung Muscle Skin Spleen Stomach Thymus Thyroid )
# $ for i in ${organs[@]}; do sbatch 12-annotate_peaks.sh $i; sleep 1s; done

source ../config.sh
organ=$1
input_parallel=6

echo $organ

# get datasets
datasets=$(awk '{print $1}' ${chrombpnet_models_keep2})
datasets_organ=( $(echo ${datasets[@]} | tr ' ' '\n' | grep ${organ}) )

datasets_to_do=$( for dataset in ${datasets_organ[@]}; do

    # check if the out files exist
    out_file="${anno_peaks_dir%/}/${dataset}__peaks_bpnet.annotated.tsv"
    if [[ ! -f "$out_file" ]]; then
        echo "$dataset"
    fi

	done | uniq )

echo "@ Processing cell types: ${datasets_to_do[@]}"

parallel --linebuffer -j ${input_parallel} bash 12-jobscript.sh {} ::: ${datasets_to_do[@]}
