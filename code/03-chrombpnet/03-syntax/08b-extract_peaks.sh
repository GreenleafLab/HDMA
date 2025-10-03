#!/usr/bin/bash


source ../config.sh

# Load required modules
module load biology samtools bedtools

datasets=$(awk '{print $1}' "${chrombpnet_models_keep2}")

region_types=( "Distal" "Promoter" "Exonic" "Intronic")

for region in ${region_types[@]}; do
    echo $region

    region_out_dir="${base_dir}/03-syntax/08/peaks/${region}"
    [[ -d "${region_out_dir}" ]] || mkdir -p "${region_out_dir}"

    # Extract peaks matching the region type
    for dataset in $datasets; do

        echo "  -> Dataset: ${dataset}"
        peaks="${anno_peaks_dir%/}/${dataset}__peaks_bpnet.annotated.tsv"

        # Check if the peaks file exists
        if [[ ! -f "${peaks}" ]]; then
          echo "Warning: Peaks file ${peaks} does not exist. Skipping..."
          continue
        fi
        
        # Create a BED file for the filtered peaks
        filtered_bed="${region_out_dir}/peaks.${region}.${dataset}.bed"

        # Check if the output file exists
        if [[ -f "${filtered_bed}" ]]; then
          echo "@ done for $motif"
          continue
        fi

        awk -F'\t' -v r="${region}" '$9 == r {print $0}' "${peaks}" | sort -k1,1 -k2,2n >> "${filtered_bed}"

    done

done

echo "#####################################################"
echo "### @ done."
echo "#####################################################"







