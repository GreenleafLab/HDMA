#!/usr/bin/bash

# Purpose: this script runs through batches of clusters of patterns (clusters
# identified by gimme cluster) and merges similar patterns

set -euo pipefail

# SETUP ------------------------------------------------------------------------

# source configuration variables
source ../config.sh

BATCH_SIZE=5
bias_params="Heart_c0_thresh0.4"
out_dir=${modisco_merged_dir}
model_head="counts"

# set inputs
cluster_key="${gimme_cluster_dir%/}/gimme_cluster_all_cluster_key.tsv"

# SUBMIT JOBS ------------------------------------------------------------------

# get batches
batches=( $(cut -f 4 $cluster_key | sort | uniq) )

for batch in ${batches[@]}; do

    job_name="03-merge_modisco_batch${batch}"
    JOBSCRIPT=03-jobscript.sh

    if [[ -f ${out_dir}/.${batch}.done ]]; then
        echo "@ found done file for batch ${batch}, skipping..."
    else

        echo "@ submitting job with: "
        echo "--out-dir ${out_dir}"
        echo "--model-head ${model_head}"
        echo "--cluster-key ${cluster_key}"
        echo "--modisco-dir ${modisco_scratch%/}/bias_${bias_params}/"
        echo "--contribs-dir ${contribs_scratch%/}/bias_${bias_params}"
        echo "--batch ${batch}"
    
        sbatch -J ${job_name} ${JOBSCRIPT} ${out_dir} ${model_head} ${cluster_key} "${modisco_scratch%/}/bias_${bias_params}" "${contribs_scratch%/}/bias_${bias_params}" ${batch}
    
        sleep 2s
    fi
done
