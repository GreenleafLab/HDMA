#!/bin/bash

# Purpose: Calculate similarity of genomic regions using the modified Jaccard index
# from bedtools, between genomic regions for some feature from each cluster.
# 
# https://bedtools.readthedocs.io/en/latest/content/tools/jaccard.html
# Favorov et al [1] reported the use of the Jaccard statistic for genome intervals:
# specifically, it measures the ratio of the number of intersecting base pairs between
# two sets to the number of base pairs in the union of the two sets. The bedtools
# jaccard tool implements this statistic, yet modifies the statistic such that the
# length of the intersection is subtracted from the length of the union. As a result,
# the final statistic ranges from 0.0 to 1.0, where 0.0 represents no overlap and 1.0
# represent complete overlap.
# # Usage:
# $ bash <script> <dir> <key>

# Check if the right number of arguments are provided
if [[ $# -ne 2 ]]; then
    echo "Usage: bash $0 bash <dir> <key>"
    exit 1
fi

source ../config.sh

module load biology bedtools

# inputs
BED_DIR=$1
KEY=$2

# list all BED files in the directory
CLUSTERS=("$BED_DIR"/*.bed)

num_clusters=${#CLUSTERS[@]}

# outputs
OUT="${base_dir}/03-syntax/08/"
RESULTS_FILE="${OUT}/jaccard_results.${KEY}.txt"

echo -e "cluster_A\tcluster_B\tlength_intersection\tlength_union\tjaccard_index\tnum_intersections" > "$RESULTS_FILE"

# loop through each pair of cluster files
# don't repeat A <--> B and B <--> A
for ((i=0; i<num_clusters-1; i++)); do
    for ((j=i; j<num_clusters; j++)); do
    
        cluster_a="${CLUSTERS[$i]}"
        cluster_b="${CLUSTERS[$j]}"

        jaccard_output=$(bedtools jaccard -a "$cluster_a" -b "$cluster_b")

        echo "$jaccard_output" | awk -v cl1="$(basename "$cluster_a")" -v cl2="$(basename "$cluster_b")" '
        {
            if (NR > 1) {
		            print cl1 "\t" cl2 "\t" $1 "\t" $2 "\t" $3 "\t" $4
	          }
        }' >> "$RESULTS_FILE"
    done
done

echo "@ done $RESULTS_FILE."
