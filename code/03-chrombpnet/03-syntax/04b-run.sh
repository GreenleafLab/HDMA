#!/usr/bin/bash

JOBSCRIPT=04b-jobscript.sh
MOTIF_FILE=04b-composites_to_test.tsv

# count motifs
NUM_JOBS=$(wc -l < "${MOTIF_FILE}")

# account for header
NUM_JOBS=$((NUM_JOBS - 1))

echo "Number of motifs to test: ${NUM_JOBS}"

# sbatch --array=1-${NUM_JOBS}%40 ${JOBSCRIPT}
sbatch --array=1-1 ${JOBSCRIPT}