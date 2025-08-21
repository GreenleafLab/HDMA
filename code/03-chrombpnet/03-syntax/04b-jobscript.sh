#!/bin/bash
#SBATCH --output=../../logs/03-chrombpnet/03/04b/%x-%A_%a.out
#SBATCH -p akundaje,wjg
#SBATCH --time=36:00:00
#SBATCH -c 2
#SBATCH --mem=40G
#SBATCH -G 1
#SBATCH --requeue
#SBATCH --open-mode=append

# set -eo pipefail

# load conda environment
eval "$(conda shell.bash hook)"
conda activate tangermeme

source ../config.sh

CLUSTER_FILE=$chrombpnet_models_keep2

# get cluster line to read
echo $SLURM_ARRAY_TASK_ID

CLUSTER=$(awk -v line="$SLURM_ARRAY_TASK_ID" 'NR==line {print $1}' "${CLUSTER_FILE}")

if [ -z "${CLUSTER}" ]; then
    echo "Error: Could not find cluster on line ${SLURM_ARRAY_TASK_ID} of ${CLUSTER_FILE}."
    exit 1
fi

echo "--- Job Details ---------------------"
echo "Job Name: ${SLURM_JOB_NAME}, Job ID: ${SLURM_JOB_ID}"
echo "Array Task ID: ${SLURM_ARRAY_TASK_ID}"
echo "Cluster: ${CLUSTER}"
echo "-------------------------------------"

echo "Executing: python -u 04b-in_silico_test_cooperativity.py --env sherlock --cluster '${CLUSTER}'"
python -u 04b-in_silico_test_cooperativity.py \
    --env sherlock \
    --cluster ${CLUSTER}
    
    # for debugging:
    #  \
    # --motif '434|SOX_SOX#1' # TODO