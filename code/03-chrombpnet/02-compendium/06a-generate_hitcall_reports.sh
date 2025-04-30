#!/bin/bash
#SBATCH --job-name=06-generate_report
#SBATCH --output=../../logs/03-chrombpnet/02/06/%x-%j.out
#SBATCH --partition=akundaje,wjg,biochem,sfgf
#SBATCH --mem-per-cpu=10G
#SBATCH --cpus-per-task=4
#SBATCH --time=04:00:00

set -euo pipefail

# Purpose: runner script to generate and execute commands for producing reports to
# evaluate FiNeMo hitcalling. Polars used in the reporting function is already
# multi-threaded so we don't further parallelize here.
# Usage: 
# $ organs=( Adrenal Brain Eye Heart Liver Lung Muscle Skin Spleen Stomach Thymus Thyroid )
# $ for i in ${organs[@]}; do sbatch 06-generate_hitcall_reports.sh $i; sleep 1s; done

source ../config.sh

organ=$1
version=0.23
alpha="0.8"
finemo_param1="counts_v${version}_a${alpha}_all"
finemo_param2="counts_v${version}_a${alpha}_nocompo"
motif_anno="04d-ChromBPNet_de_novo_motifs.tsv"
anno_drop="exclude"

# get datasets
datasets=$(awk '{print $1}' ${chrombpnet_models_keep2})

# filter to elements in the array that match the ${organ} variable
datasets_organ=( $(echo ${datasets[@]} | tr ' ' '\n' | grep ${organ}) )

datasets_to_do=$( for dataset in ${datasets_organ[@]}; do
	    
    # check that both input files exist, and the output file does not
    in_file="${hits_unified_scratch}/${dataset}/${finemo_param1}/hits.bed.gz"
    in_file2="${hits_unified_scratch}/${dataset}/${finemo_param2}/hits.bed.gz"
    out_file2="${hits_unified_scratch}/${dataset}/${finemo_param2}/custom_report/report.html"

    if [[ -f "$in_file" && -f "$in_file2" && ! -f "$out_file2" ]]; then
        echo "$dataset"
    fi
	
	done | uniq )

# DEBUG
# datasets_to_do=( Eye_c11 Liver_c4 )

echo "@ Processing cell types: ${datasets_to_do[@]}"

for dataset in ${datasets_to_do[@]}; do
   bash 06a-jobscript.sh ${dataset} ${finemo_param1} ${motif_anno} ${anno_drop} ${alpha}
   bash 06a-jobscript.sh ${dataset} ${finemo_param2} ${motif_anno} ${anno_drop} ${alpha}
done

