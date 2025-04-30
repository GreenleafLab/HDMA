#!/bin/bash
#SBATCH --job-name=06b-reconcile_hits
#SBATCH --output=../../logs/03-chrombpnet/02/06/%x-%j.out
#SBATCH --partition=akundaje,wjg,biochem,sfgf
#SBATCH --mem-per-cpu=16G
#SBATCH --cpus-per-task=6
#SBATCH --time=04:00:00

# Purpose: runner script to filter hits (remove low quality or other motifs)
# marked for exclusion, label hits with annotated motif names,
# reconcile hits between different runs and for identical
# motifs, and annotate FiNeMo hits with genomic localization information. 
#
# Usage:
# $ organs=( Adrenal Brain Eye Heart Liver Lung Muscle Skin Spleen Stomach Thymus Thyroid )
# $ for i in ${organs[@]}; do sbatch 06b-reconcile_hits.sh $i; sleep 1s; done

echo $1

set -euo pipefail

source ../config.sh

organ=$1

input_parallel=6
alpha="0.8"
finemo_param1="counts_v0.23_a${alpha}_all"
finemo_param2="counts_v0.23_a${alpha}_nocompo"
motif_annotation="04d-ChromBPNet_de_novo_motifs.tsv"
annotation_drop="exclude"

# get datasets
datasets=$(awk '{print $1}' ${chrombpnet_models_keep2})
datasets_organ=( $(echo ${datasets[@]} | tr ' ' '\n' | grep ${organ}) )

# NOTE: to update labels, first clean up previous run:

datasets_to_do=$( for dataset in ${datasets_organ[@]}; do

    # check if the in/out files exist
    # if not, keep the dataset & get list of all unique ones
    in_file="${hits_unified_scratch}/${dataset}/${finemo_param1}/seqlet_report.tsv"
    in_file2="${hits_unified_scratch}/${dataset}/${finemo_param2}/seqlet_report.tsv"
    out_file="${hits_reconciled_scratch}/${dataset}/${finemo_param1}/${dataset}__hits_unified.${finemo_param1}.reconciled.bed.gz"
    out_file2="${hits_reconciled_scratch}/${dataset}/${finemo_param1}/hits_unique.reconciled.annotated.tsv.gz"
    if [[ -f "$in_file" && -f "$in_file2" && ! -f "$out_file2" ]]; then
        echo "$dataset"
    fi

	done | uniq )

datasets_to_do=${datasets_organ[@]}

echo "@ Processing cell types: ${datasets_to_do[@]}"

for i in ${datasets_to_do[@]}; do
    echo $i
    bash 06b-jobscript.sh $i ${finemo_param1} ${finemo_param2} ${motif_annotation} ${annotation_drop}
done

parallel --linebuffer -j ${input_parallel} bash 06b-jobscript.sh {} ${finemo_param1} ${finemo_param2} ${motif_annotation} ${annotation_drop} ::: ${datasets_to_do[@]}
