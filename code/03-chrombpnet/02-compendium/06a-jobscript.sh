#!/bin/bash

set -eo pipefail

# get inputs
dataset=$1
finemo_param=$2
motif_anno=$3
anno_drop=$4
alpha=$5

# dataset="Thyroid_c7"
# finemo_param="counts_a0.7_post_filter"

# load conda environment
eval "$(conda shell.bash hook)"
conda activate finemo

source ../config.sh

# set up params
bias_params="Heart_c0_thresh0.4"
unified_h5="${base_dir%/}/02-compendium/modisco_compiled/modisco_compiled.h5"

peaks_bed="${chrombpnet_peaks_dir%/}/${dataset}__peaks_bpnet.narrowPeak"
celltype_h5="${modisco_dir%/}/bias_${bias_params}/${dataset}/counts_modisco_output.h5"
finemo_out="${hits_unified_scratch}/${dataset}/${finemo_param}/"
finemo_npz="${finemo_out%/}/intermediate_inputs.npz"

# make out dir
mkdir -p ${finemo_out}/custom_report

python 06a-unified_hitcall_report.py \
    --hits ${finemo_out}/hits.tsv \
    --regions ${finemo_npz} \
    --peaks ${peaks_bed} \
    --out-dir ${finemo_out}/custom_report/ \
    --report "06a-template.html" \
    --merged-motif-map "${modisco_comp_anno_dir}/modisco_merged_reports.tsv" \
    --merged-motif-anno ${motif_anno} \
    --modisco-region-width 400 \
    --celltype ${dataset} \
    --celltype-modisco-h5 ${celltype_h5} \
    --unified-modisco-h5 ${unified_h5} \
    --annotation-drop ${anno_drop} \
    --alpha ${alpha}