#!/bin/bash

# Purpose: Calculate the proportion of hits to positive patterns that switch to
# negative patterns in other cell types, and vice versa, using the Jaccard index.
# This would represent the proportion of base pairs that are assigned to one direction 
# of contribution in one cell type which switch to the opposite in another cell type.
#
# # Usage:
# $ bash <script> <dir_A> <dir_2> <key>

# Check if the right number of arguments are provided
if [[ $# -ne 3 ]]; then
    echo "Usage: bash $0 bash <dir_A> <dir_B> <key>"
    exit 1
fi

source ../config.sh

module load biology bedtools

# inputs
BED_DIR_A=$1
BED_DIR_B=$2
KEY=$3

# list all BED files in the directory
CLUSTERS_A=("$BED_DIR_A"/*.bed)
CLUSTERS_B=("$BED_DIR_B"/*.bed)

num_clusters=${#CLUSTERS_A[@]}

# outputs
OUT="${base_dir}/03-syntax/08/"
RESULTS_FILE="${OUT}/jaccard_results.${KEY}.txt"

echo -e "cluster_A\tcluster_B\tlength_intersection\tlength_union\tjaccard_index\tnum_intersections" > "$RESULTS_FILE"

# loop through each pair of cluster files, where A <--> B and B <--> A
# are non-commutative and both computed separately
for ((i=0; i<num_clusters-1; i++)); do
    for ((j=0; j<num_clusters; j++)); do
    
        cluster_a="${CLUSTERS_A[$i]}"
        cluster_b="${CLUSTERS_B[$j]}"

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
