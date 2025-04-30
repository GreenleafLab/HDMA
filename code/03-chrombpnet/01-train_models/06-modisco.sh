#!/bin/bash

# Purpose: this script does a deep MoDisco run with 1M seqlets for each of the
# positive and negative sets, on averaged contribution scores for counts head,
# for each cell type.

# PARAMETERS -------------------------------------------------------------------

# source configuration variables
source ../config.sh

# set bias model
bias_params="Heart_c0_thresh0.4"
bias_model="${bias_dir%/}/${bias_params}/models/bias.h5"

# set MODISCO params
max_seqlets=1000000
num_leiden=2
meme_db="${vierstra_v1_dir}/motifs.meme.txt"
num_matches=10

# set outdir on scratch for faster IO
in_dir="${contribs_scratch%/}/bias_${bias_params}"
out_dir="${modisco_scratch%/}/bias_${bias_params}"
[[ -d ${out_dir} ]] || mkdir -p ${out_dir}

JOBSCRIPT=06-jobscript.sh

# CONSTRUCT COMMANDS -----------------------------------------------------------

datasets=$(awk '{print $1}' ${chrombpnet_models_keep})

# DEBUG
# datasets=( Muscle_c7 Thymus_c16 Thyroid_c10 )

for dataset in ${datasets[@]}; do
  echo ${dataset}

  # set up common inputs/outputs
	img_suffix_dir="./"
	
	# run for COUNTS -----------------------------------------
	job_name="06-modisco_${dataset}_counts"
	counts_output_memedb=${out_dir}/${dataset}/${dataset}_memedb.counts.txt
	counts_peak_shaps=${in_dir}/${dataset}/average_shaps.counts.h5
	counts_modisco_output=${out_dir}/${dataset}/counts_modisco_output.h5
	counts_report_outdir=${out_dir}/${dataset}/counts_modisco_report
	[[ -d ${counts_report_outdir} ]] || mkdir -p ${counts_report_outdir}

  # check if outs exist already
  if [[ -f $counts_modisco_output && -f ${counts_report_outdir}/motifs.html ]]; then 
    echo "@ done ${dataset} modisco on counts."
  else
    echo "@ Running modisco with: "
  	echo -e "\t job name: ${job_name}"
  	echo -e "\t counts shaps: ${counts_peak_shaps}"
  	echo -e "\t counts modisco output: ${counts_modisco_output}"
  	echo -e "\t counts report out: ${counts_report_outdir}"
  
  	sbatch -J ${job_name} ${JOBSCRIPT} ${counts_peak_shaps} \
                                       ${max_seqlets} \
                                       ${num_leiden} \
                                       ${counts_modisco_output} \
                                       ${counts_report_outdir} \
                                       ${img_suffix_dir} \
                                       ${meme_db} \
                                       ${num_matches} \
                                       ${counts_output_memedb}

  	sleep 2s
  fi
	
done
