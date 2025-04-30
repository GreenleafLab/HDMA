#!/bin/bash

# Purpose: this script generates bigwigs of the predicted bias-corrected or uncorrected
# chromatin accessibility profiles.

# PARAMETERS -------------------------------------------------------------------

# source configuration variables
source ../config.sh

# set bias model
bias_params="Heart_c0_thresh0.4"
model_dir=${models_dir%/}/bias_${bias_params}
ref_fasta="${refs}/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta"

echo $model_dir
echo $ref_fasta
echo $chromsizes

JOBSCRIPT=08-jobscript.sh

datasets=$(awk '{print $1}' ${chrombpnet_models_keep2})

# types of predictions
modes=("bias_corrected" "uncorrected")

for dataset in "${datasets[@]}"; do
  echo "Processing dataset: ${dataset}"

  for mode in "${modes[@]}"; do
    echo "  Mode: ${mode}"

    # construct paths based on the current mode
    if [[ "$mode" == "bias_corrected" ]]; then
      out_dir="${preds_scratch%/}/bias_corrected"
      final_out_dir="${preds_dir%/}/bias_corrected"
      out_key="nobias"
      model_suffix="_nobias"
    else
      out_dir="${preds_scratch%/}/uncorrected"
      final_out_dir="${preds_dir%/}/uncorrected"
      out_key="uncorrected"
      model_suffix=""
    fi

    # specify models
    model_fold_0="${model_dir}/${dataset}/fold_0/models/chrombpnet${model_suffix}.h5"
    model_fold_1="${model_dir}/${dataset}/fold_1/models/chrombpnet${model_suffix}.h5"
    model_fold_2="${model_dir}/${dataset}/fold_2/models/chrombpnet${model_suffix}.h5"
    model_fold_3="${model_dir}/${dataset}/fold_3/models/chrombpnet${model_suffix}.h5"
    model_fold_4="${model_dir}/${dataset}/fold_4/models/chrombpnet${model_suffix}.h5"

    peaks_file="${chrombpnet_peaks_dir%/}/${dataset}__peaks_bpnet.narrowPeak"
    out_prefix="${out_dir%/}/${dataset}_avg"
    final_out_prefix="${final_out_dir%/}/${dataset}_avg"
    final_out_file1="${final_out_prefix}_chrombpnet_${out_key}_preds_w_logcounts.bed"
    final_out_file2="${final_out_prefix}_chrombpnet_${out_key}.bw"

    if [[ -f "${final_out_file1}" && -f "${final_out_file2}" ]]; then
      echo -e "\t\tfound completed ${mode} predictions for ${dataset}, skipping..."
    else
      echo "Generating predictions for ${dataset} (${mode})"
      echo "Final Output File 1: ${final_out_file1}"
      echo "Final Output File 2: ${final_out_file2}"

      job_name="08-predict_${dataset}_${mode}"
      echo "@ generating ${out_prefix} ${out_key} predictions"

      echo "sbatch -J ${job_name} ${JOBSCRIPT} ${dataset}"
      echo " ${peaks_file}"
      echo " ${ref_fasta}"
      echo " ${chromsizes}"
      echo " ${out_prefix}"
      echo " ${out_key}"
      echo " ${model_fold_0}"
      echo " ${model_fold_1}"
      echo " ${model_fold_2}"
      echo " ${model_fold_3}"
      echo " ${model_fold_4}"

      sbatch -J "${job_name}" "${JOBSCRIPT}" "${dataset}" \
          "${peaks_file}" \
          "${ref_fasta}" \
          "${chromsizes}" \
          "${out_prefix}" \
          "${out_key}" \
          "${model_fold_0}" \
          "${model_fold_1}" \
          "${model_fold_2}" \
          "${model_fold_3}" \
          "${model_fold_4}"

      sleep 3s
    fi
  done
done

