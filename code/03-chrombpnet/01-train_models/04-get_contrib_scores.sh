#!/bin/bash

# Purpose: this script calls the chrombpnet helper to get contribution 
# score bigwigs, as computed with DeepLIFT: https://github.com/kundajelab/chrombpnet/wiki/Generate-contribution-score-bigwigs
# Based off code from Ryan Zhao / Salil Deshpande.

# PARAMETERS -------------------------------------------------------------------

# source configuration variables
source ../config.sh

# set bias model
bias_params="Heart_c0_thresh0.4"
bias_model="${bias_dir%/}/${bias_params}/models/bias.h5"

# set parameters
ref_fasta="${refs}/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta"
out_dir="${contribs_scratch%/}/bias_${bias_params}/"
[[ -d ${out_dir} ]] || mkdir -p ${out_dir}



# CONSTRUCT COMMANDS -----------------------------------------------------------

# for every [cluster]__sorted.tsv file, extract [cluster]
datasets=$(awk '{print $1}' ${chrombpnet_models_keep})

# DEBUG
# datasets=( Muscle_c7 Thymus_c16 Thyroid_c10 )

for dataset in ${datasets[@]}; do
  echo ${dataset}

	peaks_file=${chrombpnet_peaks_dir%/}/${dataset}__peaks_bpnet.narrowPeak
	dataset_dir=${out_dir%/}/${dataset}
	[[ -d ${dataset_dir} ]] || mkdir ${dataset_dir}

	for fold in {0..4}; do
		fold_name=fold_${fold}
		echo -e "\t${fold_name}"

		model_file=${models_dir%/}/bias_${bias_params}/${dataset}/${fold_name}/models/chrombpnet_nobias.h5
		
		[[ -f ${model_file} ]] || echo -e "\t\tModel file missing"
		
		fold_dir=${dataset_dir}/${fold_name}

		[[ -d ${fold_dir} ]] || mkdir ${fold_dir}

		out_prefix=${fold_dir}/peaks_shap

		if [[ -f "${out_prefix}.counts_scores.bw" && -f "${out_prefix}.profile_scores.bw" ]]; then
			echo -e "\t\tfound completed shap scores, skipping..."
		else 
			job_name="04-get_contrib_${dataset}_${fold_name}"
			JOBSCRIPT=04-jobscript.sh
			
			echo "Running chrombpnet contribs_bw"
			echo "--genome ${ref_fasta}"
			echo "--regions ${peaks_file}"
			echo "--model-h5 ${model_file}"
			echo "--output-prefix ${out_prefix}"
			echo "--chrom-sizes ${chrom_sizes}"
			
			sbatch -J ${job_name} ${JOBSCRIPT} ${ref_fasta} \
				${peaks_file} \
				${model_file} \
				${out_prefix} \
				${chromsizes}
			
			sleep 5s
		fi
	done
	sleep 20s
done
