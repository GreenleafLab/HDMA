#!/usr/bin/bash

# source configuration variables
source ../config.sh

# count peaks
wc -l $chrombpnet_peaks_dir/*.narrowPeak > ${base_dir%/}/01-models/qc/npeaks.tsv
