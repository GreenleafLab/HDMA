#!/bin/bash


# Purpose: for all the cell types in each organ, run gimmemotifs cluster routine
# for all motifs from those clusters. The output is a set of PFMs, one per cluster, which are the averages
# of the PFMs for all motifs in the respective cluster.

# NOTE: re: gimme motifs installation
# I installed gimmemotifs, inside a _conda_ environment but using mamba install.
# mamba install let me successfully install the package, while the conda env
# lets me activate it / access gimme inside a shell as below.
# This is run for both pos_patterns and neg_patterns separately.
	
set -euo pipefail

# SETUP ------------------------------------------------------------------------

# source configuration variables
source ../config.sh

# RUN GIMME CLUSTER ------------------------------------------------------------

# for each organ, run gimme cluster on the positive patterns

prefixes=( Adrenal Brain Eye Heart Liver Lung Muscle Skin Spleen Stomach Thymus Thyroid )

for prefix in ${prefixes[@]}; do

  key="${prefix}.counts.neg_patterns"
	t=0.8
	ncpus="16"

	params="t${t}_n${ncpus}"
	out_dir=${gimme_cluster_dir%/}/${key}_${params}

	echo "@ input: ${pfm_dir%/}/${key}.pfm"
	echo "@ output: ${out_dir}"

	# make output directory if it doesn't exist
	[[ -d ${out_dir} ]] || mkdir -p ${out_dir}

	# RUN GIMME CLUSTER ------------------------------------------------------------

	job_name="02-gimme_cluster_${key}_${params}"
	JOBSCRIPT=02-jobscript.sh
	echo "@ submitting: gimme cluster ${pfm_dir%/}/${key}.pfm ${out_dir} ${t}"

	sbatch -J ${job_name} ${JOBSCRIPT} ${pfm_dir%/}/${key}.pfm ${out_dir} ${t}

	sleep 2s

done
