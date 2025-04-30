#!/usr/bin/bash
#SBATCH --job-name="06-make_bw"
#SBATCH --time=12:00:00
#SBATCH --output=../../logs/03-chrombpnet/00/%x-%j.out
#SBATCH --partition=akundaje,owners
#SBATCH --cpus-per-task=6
#SBATCH --mem=30G

# Create unnormalized bigwigs from the observed pseudobulk accessibility per cluster.

# load conda environment
eval "$(conda shell.bash hook)"
conda activate chrombpnet

# source configuration variables
source ../config.sh

input_parallel=6
export out_dir="${bigwigs_scratch%/}"
export script_loc="${chrombpnet_code}/chrombpnet/helpers/preprocessing/reads_to_bigwig.py"
export fasta_file="${refs}/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta"
export cluster_frags_dir="${cluster_frags_dir%/}"
export chromsizes=$chromsizes

makebigwig () {
  dataset=${1}
	frag_file="${cluster_frags_dir}/fragments/${dataset}__sorted.tsv"
	out_file="${out_dir}/${dataset}"
	
	head $frag_file
	echo $out_file
	
	echo "${dataset} starting"
	python3.8 ${script_loc} -g ${fasta_file} -ifrag ${frag_file} -c ${chromsizes} -o ${out_file} -d ATAC
	echo "${dataset} done"
}
export -f makebigwig

all_pseudorep_dir="${cluster_frags_dir%/}/fragments"

# for every [cluster]__sorted.tsv file, extract [cluster]
all_sorted=$(ls $all_pseudorep_dir/*__sorted.tsv | xargs -n 1 -I {} basename {} __sorted.tsv)	

datasets=$( for dataset in ${all_sorted[@]}; do
	    
    # check if the done file exists
    # if not, keep the dataset & get list of all unique ones
    done_file="${bigwigs_dir%/}/${dataset}_unstranded.bw"
    if [[ ! -f "$done_file" ]]; then
        echo "$dataset"
    fi
	
	done | uniq )


echo "@ Producing bigwigs for: ${datasets}"

# DEBUG:
# makebigwig Eye_c0

parallel -j ${input_parallel} makebigwig {} ::: ${datasets}

echo "@ done"
