#!/bin/bash
#SBATCH --output=../../logs/03-chrombpnet/02/02/02b-gimme_cluster_all_%j.out
#SBATCH -p akundaje,wjg,sfgf,biochem
#SBATCH -t 02:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=60G

# Purpose: use gimme to cluster the clusters averages previously identified in Step 02a.

set -euo pipefail

# source configuration variables
source ../config.sh

# load conda env	
eval "$(conda shell.bash hook)"
conda activate gimme

# ------------------------------------------------------------------------------
# POS PATTERNS
# first, concatenate the gimme clusters from each cell type
organs=( Adrenal Brain Eye Heart Liver Lung Muscle Skin Spleen Stomach Thymus Thyroid )

# clear the file if it exists already
all_pfm="${pfm_dir%/}/all_clustered_clusters.counts.neg_patterns.pfm"
[[ -f $all_pfm ]] && rm $all_pfm

for organ in ${organs[@]}; do

	key="${organ}.counts.neg_patterns"
	t=0.8
	ncpus="16"
	params="t${t}_n${ncpus}"
	gimme_out=${gimme_cluster_dir%/}/${key}_${params}/clustered_motifs.pfm
	echo ${gimme_out}
	# prefix the cluster average motif names with the organ name, then concatenate
	if [[ -f $gimme_out ]]; then
		sed "s/>/>${organ}__/g" $gimme_out >> ${all_pfm}
	fi

done

# number of motifs
grep -E '>' ${all_pfm} | wc -l

# ------------------------------------------------------------------------------
# now run gimme on that new concatenated set of PFMs
key="all_clustered_clusters.counts.neg_patterns"
t=0.8
ncpus="16"
params="t${t}_n${ncpus}"

out_dir=${gimme_cluster_dir%/}/${key}_${params}
input=${all_pfm}

[[ ! -d $out_dir ]] && mkdir -p $out_dir

echo "@ input: ${input}"
echo "@ output: ${out_dir}"

echo "@ running: gimme cluster ${input} ${out_dir} ${t}"

gimme cluster ${input} ${out_dir} -t ${t} -N 16

echo "@ done."


# ------------------------------------------------------------------------------
# NEG PATTERNS
# first, concatenate the gimme clusters from each cell type (only the organs that have negative patterns)
organs=( Brain Eye Heart Lung Muscle Skin Spleen Thymus )

# clear the file if it exists already
all_pfm="${pfm_dir%/}/all_clustered_clusters.counts.neg_patterns.pfm"
[[ -f $all_pfm ]] && rm $all_pfm

for organ in ${organs[@]}; do

	key="${organ}.counts.neg_patterns"
	t=0.8
	ncpus="16"
	params="t${t}_n${ncpus}"
	gimme_out=${gimme_cluster_dir%/}/${key}_${params}/clustered_motifs.pfm
	echo ${gimme_out}
	# prefix the cluster average motif names with the organ name, then concatenate
	if [[ -f $gimme_out ]]; then
		sed "s/>/>${organ}__/g" $gimme_out >> ${all_pfm}
	fi

done

# Liver, Stomach only had one neg pattern across cell types in the organ; concatenate it here
sed "s/>/>Liver__/g" ${pfm_dir%/}/Liver.counts.neg_patterns.pfm >> ${all_pfm}
sed "s/>/>Stomach__/g" ${pfm_dir%/}/Stomach.counts.neg_patterns.pfm >> ${all_pfm}

# number of motifs
grep -E '>' ${all_pfm} | wc -l

# ------------------------------------------------------------------------------
# now run gimme on that new concatenated set of PFMs
key="all_clustered_clusters.counts.neg_patterns"
t=0.8
ncpus="16"
params="t${t}_n${ncpus}"
	
out_dir=${gimme_cluster_dir%/}/${key}_${params}
input=${all_pfm}

[[ ! -d $out_dir ]] && mkdir -p $out_dir

echo "@ input: ${input}"
echo "@ output: ${out_dir}"

echo "@ running: gimme cluster ${input} ${out_dir} ${t}"

gimme cluster ${input} ${out_dir} -t ${t} -N 16

echo "@ done."
