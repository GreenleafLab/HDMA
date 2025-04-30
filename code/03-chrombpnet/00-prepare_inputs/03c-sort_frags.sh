#!/usr/bin/bash

# Purpose: sort the fragemnts for each cluster.
# Adapted from Kundaje Lab / Salil Deshpande. This is run on Kundaje lab cluster.

# source configuration variables
source ../config.sh

# specify inuts
input_basedir="${base_dir%/}/00-inputs/"
frags_dir=$cluster_frags_dir
input_parallel=6
input_tmp=$labcluster_scratch

echo $input_basedir
echo $frags_dir

# begin script
# - nice processes to avoid hogging cluster
# - explicitly limit memory usage by specifying buffersize using -S
# - use an alternate temp directory in /srv/scratch to avoid out-of-space issues with /tmp
dosort () {
	dataset=${1}
	frags_dir=${2}
	input_tmp=${3}

	echo "${dataset} sorting fragments..."
	nice -n 19 sort -k 1,1 -k 2,2n -S 20G -T ${input_tmp} ${frags_dir}/fragments/${dataset}.tsv -o ${frags_dir}/fragments/${dataset}__sorted.tsv
	echo "${dataset} sorting pseudorepT..."
	nice -n 19 sort -k 1,1 -k 2,2n -S 20G -T ${input_tmp} ${frags_dir}/pseudorepT/${dataset}.tsv -o ${frags_dir}/pseudorepT/${dataset}__sorted.tsv
	echo "${dataset} sorting pseudorep1..."
	nice -n 19 sort -k 1,1 -k 2,2n -S 20G -T ${input_tmp} ${frags_dir}/pseudorep1/${dataset}.tsv -o ${frags_dir}/pseudorep1/${dataset}__sorted.tsv
	echo "${dataset} sorting pseudorep2..."
	nice -n 19 sort -k 1,1 -k 2,2n -S 20G -T ${input_tmp} ${frags_dir}/pseudorep2/${dataset}.tsv -o ${frags_dir}/pseudorep2/${dataset}__sorted.tsv

	touch "${frags_dir}/fragments/.${dataset}.done"

	echo "@@ done ${dataset}"
}
export -f dosort

all_fragments_dir="${frags_dir%/}/fragments"
all_fragments=$(ls ${all_fragments_dir} | xargs -n 1 -I {} basename {} __sorted.tsv)

datasets=$( for frag_file in ${all_fragments[@]}; do

    # extract [cluster] from filename
    dataset=${frag_file%.tsv}

    # check if the done file exists; if not, keep the dataset & get list of all unique ones
    done_file="${frags_dir}/fragments/.${dataset}.done"
    if [[ ! -f "$done_file" ]]; then
        echo "$dataset"
    fi

	done | uniq )


echo "@ running sort on clusters: ${datasets}"

# for dataset in ${datasets[@]}; do dosort $dataset $frags_dir; done
parallel --linebuffer -j ${input_parallel} dosort {} $frags_dir $input_tmp ::: ${datasets}

echo "@ done"
