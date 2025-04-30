#!/usr/bin/bash
#SBATCH --job-name="03a-cat"
#SBATCH --time=04:00:00
#SBATCH --output=../../logs/03-chrombpnet/00/%x-%j.out
#SBATCH --partition=akundaje,owners
#SBATCH --cpus-per-task=6
#SBATCH --mem=40G

# Purpose: concatenate the fragments across samples for each cluster.
# Fragments are then sorted in Step 3b.
# Adapted from Kundaje Lab / Salil Deshpande.

# source configuration variables
source ../config.sh

# specify inuts
input_basedir="${base_dir%/}/00-inputs/"
input_parallel=6

in_dir=$sample_frags_dir
out_dir=$cluster_frags_dir

echo $input_basedir
echo ${in_dir}
echo ${out_dir}

# begin script
docat () {
	dataset=${1}
	in_dir=${2}
	out_dir=${3}

	# DEBUG:
	# echo ${dataset}
	# echo ${in_dir}
	# echo ${out_dir}

	echo "${dataset} concatenating across samples..."
	find "${in_dir}/fragments" -name "${dataset}*.tsv" -exec cat {} + > "${out_dir}/fragments/${dataset}.tsv"
	find "${in_dir}/pseudorepT" -name "${dataset}*.tsv" -exec cat {} + > "${out_dir}/pseudorepT/${dataset}.tsv"
	find "${in_dir}/pseudorep1" -name "${dataset}*.tsv" -exec cat {} + > "${out_dir}/pseudorep1/${dataset}.tsv"
	find "${in_dir}/pseudorep2" -name "${dataset}*.tsv" -exec cat {} + > "${out_dir}/pseudorep2/${dataset}.tsv"

	echo "@@ done ${dataset}"
}
export -f docat

# get the unique datasets.
# NOTE: "dataset" here refers to a cluster

# for every [cluster]__[sample].tsv file,
all_fragments_dir="${in_dir%/}/fragments"
datasets=$( for frag_file in $(ls "${all_fragments_dir}"); do
	
    # extract [cluster] from filename
    dataset=${frag_file%__*.tsv}
    
    # check if the done file exists
    # if not, keep the dataset & get list of all unique ones
    done_file="${out_dir}/fragments/.${dataset}.done"
    if [[ ! -f "$done_file" ]]; then
        echo "$dataset"
    fi
	
	done | uniq )


echo "@ running cat on clusters: ${datasets}"
parallel --linebuffer -j ${input_parallel} docat {} ${in_dir} ${out_dir} ::: ${datasets}

# DEBUG:
# docat Brain_c12 ${in_dir} ${out_dir}

echo "@ done"
