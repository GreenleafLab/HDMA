#!/bin/bash

# Purpose: This script uses the finemo_gpu package by Austin Wang
# (https://github.com/austintwang/finemo_gpu) to call motif hits (referred to
# in the manuscript as predictive instances) per cell type.
#
# In this case, we use the contribution scores derived by DeepLIFT, on peaks called
# for each cell type and used for model training.
#
# However, we use the unified motif set, so that the same set
# of patterns are used to call hits in every cell type. Hit calling is run once
# without composite motifs, and once with composite motifs, thus this step depends
# to some degree on the annotations, but only for desginating composites vs
# non-composites (not for the annotation names themselves).

set -euo pipefail

# SETUP ------------------------------------------------------------------------

# source configuration variables
source ../config.sh

output_dir=${hits_unified_scratch}
mkdir -p ${output_dir}

bias_params="Heart_c0_thresh0.4"
organ="Spleen"

# params for hit calling
modisco_h5="${base_dir%/}/02-compendium/modisco_compiled/modisco_compiled.h5"
finemo_version=0.23
alpha=0.8
finemo_params="counts_v${finemo_version}_a${alpha}_all"

alpha_nocompo=0.8
finemo_params_nocompo="counts_v${finemo_version}_a${alpha_nocompo}_nocompo"
motifs_nocompo="${modisco_comp_anno_dir%/}/motif_names.non_composites.tsv"

JOBSCRIPT=05-jobscript.sh


# CONSTRUCT COMMANDS -----------------------------------------------------------

# find which cell types to keep
datasets=$(awk '{print $1}' ${chrombpnet_models_keep2})
  
# filter to elements in the array that match the ${organ} variable
datasets_organ=( $(echo ${datasets[@]} | tr ' ' '\n' | grep ${organ}) )

for dataset in ${datasets_organ[@]}; do

    echo "@ ${dataset}"

    peaks_bed=${chrombpnet_peaks_dir%/}/${dataset}__peaks_bpnet.narrowPeak

    counts_shaps_h5=${contribs_dir%/}/bias_${bias_params}/${dataset}/average_shaps.counts.h5
    counts_finemo_out="${output_dir%/}/${dataset}/${finemo_params}/"
    counts_finemo_out_nocompo="${output_dir}/${dataset}/${finemo_params_nocompo}/"

    # generate the output dirs if it doesn't exist yet
    [[ -f $counts_finemo_out ]] || mkdir -p ${counts_finemo_out}
    [[ -f $counts_finemo_out_nocompo ]] || mkdir -p ${counts_finemo_out_nocompo}
    
    if [[ -f "${counts_finemo_out%/}/hits.bed.gz" && -f "${counts_finemo_out_nocompo%/}/hits.bed.gz" ]]; then
      echo -e "\tfound hit calls, skipping..."
    else 
        job_name="05-hits_${dataset}"

        echo -e "\t@ peaks ${peaks_bed}"
        echo -e "\t@ modisco h5 ${modisco_h5}"
        echo -e "\t@ shaps h5 ${counts_shaps_h5}"
        echo -e "\t@ out1 ${counts_finemo_out%/}/hits.bed.gz"
        echo -e "\t@ out2 ${counts_finemo_out_nocompo%/}/hits.bed.gz"
        echo -e "\t@ alpha ${alpha}"
        echo -e "\t@ alpha_nocompo ${alpha_nocompo}"
        echo -e "\t@ motifs_nocompo ${motifs_nocompo}"
        
        sbatch -J ${job_name} ${JOBSCRIPT} ${dataset} \
                            ${peaks_bed} \
                            ${modisco_h5} \
                            ${counts_shaps_h5} \
                            ${counts_finemo_out} \
                            ${counts_finemo_out_nocompo} \
                            counts \
                            ${alpha} \
                            ${alpha_nocompo} \
                            ${motifs_nocompo}
    fi

    sleep 2s

done
