#!/bin/bash

set -eo pipefail

# get inputs
dataset=$1

source ../config.sh
export R_LIBS_USER=$RENV_SJ

peaks_file=${chrombpnet_peaks_dir%/}/${dataset}__peaks_bpnet.narrowPeak
out_tsv="${anno_peaks_dir%/}/${dataset}__peaks_bpnet.annotated.tsv"

echo "@ writing reconciled hits to ${out_tsv}"

Rscript 12-annotate_peaks.R ${peaks_file} ${out_tsv}