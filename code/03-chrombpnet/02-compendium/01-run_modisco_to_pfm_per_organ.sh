#!/bin/bash

# Purpose: convert the modisco CWMs to PFMs for gimme cluster. Patterns are 
# concatenated for all cell types in each organ.

set -euo pipefail

# PARAMETERS -------------------------------------------------------------------

# source configuration variables
source ../config.sh

# load conda env	
eval "$(conda shell.bash hook)"
conda activate chrombpnet

# set bias model
bias_params="Heart_c0_thresh0.4"
bias_model="${bias_dir%/}/${bias_params}/models/bias.h5"

# BEGIN SCRIPT -----------------------------------------------------------------
all_pseudorep_dir="${cluster_frags_dir%/}/fragments"

organs=( Adrenal Brain Eye Heart Liver Lung Muscle Skin Spleen Stomach Thymus Thyroid )

for organ in ${organs[@]}; do
  
  echo $organ

  key1="${organ}.counts.pos_patterns"
  out_file1=${pfm_dir%/}/${key1}.pfm
  echo $out_file1

  key2="${organ}.counts.neg_patterns"
  out_file2=${pfm_dir%/}/${key2}.pfm
  echo $out_file2

  # remove the old config file if it exists
  config=${pfm_dir%/}/${organ}.counts.config.tsv
  [[ -f $config ]] && rm $config
  
  # find which cell types to keep
  datasets=$(awk '{print $1}' ${chrombpnet_models_keep2})
  
  # filter to elements in the array that match the ${organ} variable
  datasets_organ=( $(echo ${datasets[@]} | tr ' ' '\n' | grep ${organ}) )

  for dataset in ${datasets_organ[@]}; do
  
    echo ${dataset}
    modisco_counts_h5=${modisco_scratch%/}/bias_${bias_params}/${dataset}/counts_modisco_output.h5
 
    # indicate the modisco file to use 
    if [[ ! -f $modisco_counts_h5 ]]; then
      echo "@ missing ${dataset} COUNTS modisco output."
    else 
      echo -e "${dataset}\t${modisco_counts_h5}" >> $config
    fi

  done

  # run the modisco to pfm conversion
  python 01-modisco_to_pfm.py -c $config -o $out_file1 -p pos_patterns
  python 01-modisco_to_pfm.py -c $config -o $out_file2 -p neg_patterns

done
