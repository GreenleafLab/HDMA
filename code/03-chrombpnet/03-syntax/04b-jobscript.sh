#!/bin/bash
#SBATCH --output=../../logs/03-chrombpnet/03/04b/%x-%A_%a.out
#SBATCH -p akundaje,wjg,owners,gpu
#SBATCH --time=24:00:00
#SBATCH -c 1
#SBATCH --mem=40G
#SBATCH -G 1
#SBATCH --requeue
#SBATCH --open-mode=append

# set -eo pipefail

# load conda environment
eval "$(conda shell.bash hook)"
conda activate tangermeme

source ../config.sh

MOTIF_FILE=04b-composites_to_test.tsv

# get motif line to read (add 1 to skip header)
echo $SLURM_ARRAY_TASK_ID
motif_line=$((SLURM_ARRAY_TASK_ID + 1))

motif=$(awk -v line="$motif_line" 'NR==line {print $1}' "${MOTIF_FILE}")

if [ -z "${motif}" ]; then
    echo "Error: Could not find motif on line ${motif_line} of ${MOTIF_FILE}."
    exit 1
fi

echo "--- Job Details ---------------------"
echo "Job Name: ${SLURM_JOB_NAME}, Job ID: ${SLURM_JOB_ID}"
echo "Array Task ID: ${SLURM_ARRAY_TASK_ID}"
echo "Motif: ${motif}"
echo "-------------------------------------"

echo "Executing: python -u 04b-in_silico_test_cooperativity.py --env sherlock --motif '${motif}'"
python -u 04b-in_silico_test_cooperativity.py \
    --env sherlock \
    --motif ${motif}