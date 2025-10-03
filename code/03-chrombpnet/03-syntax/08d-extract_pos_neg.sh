#!/usr/bin/bash

source ../config.sh

# Load required modules
module load biology samtools bedtools

finemo_param="counts_v0.23_a0.8_all"

datasets=$(awk '{print $1}' "${chrombpnet_models_keep2}")

# Loop over the two main pattern types you are processing
for pattern_type in "pos_patterns" "neg_patterns"; do

    echo "#####################################################"
    echo "### @ ${pattern_type}"
    echo "#####################################################"

    echo
    echo "--- @ ${motif} for ${pattern_type} ---"

    # Create a specific output directory for the current motif
    motif_out_dir="${base_dir}/03-syntax/08/patterns/${pattern_type}/"
    mkdir -p "${motif_out_dir}"

    echo $motif_out_dir

    # Extract hits for the current motif
    for dataset in $datasets; do
      echo "  -> Dataset: ${dataset}"
      hits_bed="${hits_reconciled_dir}/${dataset}/${finemo_param}/${dataset}__hits_unified.${finemo_param}.reconciled.bed.gz"
      
      # Create a BED file for the specific hits for the current motif
      cluster_bed="${motif_out_dir}/patterns.${pattern_type}.${dataset}.bed"
  
      # Check if the output file exists
      if [[ -f "${cluster_bed}" ]]; then
        echo "@ done for $motif"
        continue
      fi

      # Check if the hits file exists
      if [[ ! -f "${hits_bed}" ]]; then
        echo "Warning: Hits file ${hits_bed} does not exist. Skipping..."
        continue
      fi
      zcat "${hits_bed}" | awk -v pt="${pattern_type}" '$0 ~ pt {print $0}' | sort -k1,1 -k2,2n >> "${cluster_bed}"

    done
done

echo "#####################################################"
echo "### @ done."
echo "#####################################################"
