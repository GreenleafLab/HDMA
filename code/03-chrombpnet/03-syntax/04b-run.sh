#!/usr/bin/bash

source ../config.sh

JOBSCRIPT=04b-jobscript.sh
CLUSTER_FILE=$chrombpnet_models_keep2

# count clusters
NUM_JOBS=$(wc -l < "${CLUSTER_FILE}")

echo "Number of clusters to test: ${NUM_JOBS}"

sbatch --array=1-${NUM_JOBS}%50 ${JOBSCRIPT}
# sbatch --array=163-163 ${JOBSCRIPT}
