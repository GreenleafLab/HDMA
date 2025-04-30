#!/bin/bash

# Purpose: this script generates bigwigs based on the average contribution scores
# across folds, for each cell type.

# PARAMETERS -------------------------------------------------------------------

# source configuration variables
source ../config.sh

# set bias model
bias_params="Heart_c0_thresh0.4"
bias_model="${bias_dir%/}/${bias_params}/models/bias.h5"

# set outdir on scratch for faster IO
in_dir="${contribs_scratch%/}/bias_${bias_params}"

JOBSCRIPT=07-jobscript.sh

# CONSTRUCT COMMANDS -----------------------------------------------------------

# for every [cluster]__sorted.tsv file, extract [cluster]
datasets=$(awk '{print $1}' ${chrombpnet_models_keep2})

# DEBUG
# datasets=( Muscle_c7 Thymus_c16 Thyroid_c10 )

for dataset in ${datasets[@]}; do

	echo ${dataset}
	
	peaks_file=${chrombpnet_peaks_dir%/}/${dataset}__peaks_bpnet.narrowPeak
  
	counts_peak_shaps=${in_dir%/}/${dataset}/average_shaps.counts.h5
	counts_out_prefix=${in_dir%/}/${dataset}/average_shaps.counts
	
	profile_peak_shaps=${in_dir%/}/${dataset}/average_shaps.profile.h5
	profile_out_prefix=${in_dir%/}/${dataset}/average_shaps.profile
	
	if [[ -f "${counts_out_prefix}.bw" && -f "${profile_out_prefix}.bw" ]]; then
	  echo -e "@ Found bigwigs for ${dataset}; skipping."
  else 
  	
  	job_name="07-contrib_bw_${dataset}"
  	echo "@ Generating bigwigs for: "
  	echo -e "\t counts: ${counts_peak_shaps} to ${counts_out_prefix}.bw"
  	echo -e "\t profile: ${profile_peak_shaps} to ${profile_out_prefix}.bw"
  	echo -e "\t using ${peaks_file}"
  	sbatch -J ${job_name} ${JOBSCRIPT} ${dataset} \
                                       ${counts_peak_shaps} \
                                       ${counts_out_prefix} \
                                       ${profile_peak_shaps} \
                                       ${profile_out_prefix} \
                                       ${peaks_file} \
                                       ${chromsizes}
  fi
  
  sleep 5s

done

