#!/bin/bash

# load conda environment
eval "$(conda shell.bash hook)"
conda activate chrombpnet

ref_fasta=${1}
chrom_sizes=${2}
blacklist=${3}
peak_file=${4}
input_len=${5}
stride=${6}
fold_out=${7}
split_file=${8}

chrombpnet prep nonpeaks --genome ${ref_fasta} \
                         --chrom-sizes ${chrom_sizes} \
                         --blacklist-regions ${blacklist} \
                         --peaks ${peak_file} \
                         --inputlen ${input_len} \
                         --stride ${stride} \
                         --output-prefix ${fold_out} \
                         --chr-fold-path ${split_file}
                                                  
                                                  