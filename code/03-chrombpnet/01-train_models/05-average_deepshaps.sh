#!/bin/bash

# Purpose: average the Shap scores across folds, for counts and profile
# contribution scores separately.

# PARAMETERS -------------------------------------------------------------------

# source configuration variables
source ../config.sh

# set bias model
bias_params="Heart_c0_thresh0.4"
bias_model="${bias_dir%/}/${bias_params}/models/bias.h5"
out_dir="${contribs_scratch%/}/bias_${bias_params}/"

echo "@ contribs dir: "
echo ${out_dir}

# CONSTRUCT COMMANDS -----------------------------------------------------------

# for every line in chrombpnet_models_keep.tsv, extract the model name and folds to keep
datasets=$(awk /'{print $1}' ${chrombpnet_models_keep})

# DEBUG
# datasets=( Muscle_c7 Thymus_c16 Thyroid_c10 )

# organ="Heart"
# datasets_organ=( $(echo ${datasets[@]} | tr ' ' '\n' | grep ${organ}) )

for dataset in ${datasets[@]}; do
  
  echo ${dataset}
  dataset_dir=${out_dir%/}/${dataset}
  echo ${dataset_dir}

  # grep the models_keep file for folds to keep
  folds_keep=`grep -E "${dataset}\s" ${chrombpnet_models_keep} | awk '{print $2}'`

  echo ${folds_keep}
  
  job_name="05-avg_${dataset}"
  JOBSCRIPT=05-jobscript.sh
 
  if [[ -f "${dataset_dir}/average_shaps.counts.h5" && -f "${dataset_dir}/average_shaps.profile.h5" ]]; then
    echo "@ done ${dataset} averaging."
  else 
    echo "Submitting python 05-average_deepshaps.py ${dataset_dir} ${folds_keep} {counts,profile}"
    sbatch -J ${job_name} ${JOBSCRIPT} ${dataset_dir} ${folds_keep}
    sleep 5s
  fi

done
