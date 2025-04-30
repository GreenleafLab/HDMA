#!/bin/bash
#SBATCH --job-name=09-run_nucleoatac
#SBATCH --output=../../logs/03-chrombpnet/02/09/%x-%j.out
#SBATCH --partition=akundaje,wjg,biochem,sfgf
#SBATCH --mem-per-cpu=20G
#SBATCH --cpus-per-task=1
#SBATCH --time=06:00:00

# Purpose: This script runs NucleoATAC (https://github.com/sjessa/NucleoATAC),
# package to infer nucleosome positions from ATACseq data. This was authored
# by Alicia Schep in the Greenleaf but adapted to take fragments files as input
# instead of BAM files.
#
# Uses NucleoATAC v0.4.1:
# https://github.com/sjessa/NucleoATAC/commit/3076551c3ad26c9251c6c309d93a184e6901c600
#
# After running NucleoATAC, we also calculate the distribution of pairwise 
# distributions between nucleosome dyads, and the distribution of distances
# between motif hits and nucleosome dyads.

# set -euo pipefail

# SETUP ------------------------------------------------------------------------
# source configuration variables
source ../config.sh

finemo_param="counts_v0.23_a0.8_all"
ref_fasta="${refs}/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta"


# make out dir
output_dir=${nucleoatac_dir}
mkdir -p ${nucleoatac_dir}

JOBSCRIPT=07-jobscript.sh


# CONSTRUCT COMMANDS -----------------------------------------------------------

# find which cell types to keep
datasets=$(awk '{print $1}' ${chrombpnet_models_keep2})
  

for dataset in ${datasets[@]}; do

    echo "@ ${dataset}"

    frag_file=${cluster_frags_dir%/}/fragments/${dataset}__sorted.tsv.gz
    peaks_bed=${chrombpnet_peaks_dir%/}/${dataset}__peaks_bpnet.narrowPeak
    nucleoatac_out="${output_dir%/}/${dataset}/"
    nucleoatac_out_prefix=${nucleoatac_out}/${dataset}
    hits="${hits_reconciled_dir}/${dataset}/${finemo_param}/hits_unique.reconciled.annotated.tsv.gz"

    # generate the output dirs if it doesn't exist yet
    [[ -f ${nucleoatac_out} ]] || mkdir -p ${nucleoatac_out}
    
    if [[ -f "${nucleoatac_out%/}/${dataset}.nucmap_combined.bed.gz" ]]; then
      echo -e "\tfound nucleosome calls, skipping..."
    else
        job_name="07-nucleoatac_${dataset}"

        echo -e "\t@ dataset ${dataset}"
        echo -e "\t@ fragments ${frag_file}"
        echo -e "\t@ peaks ${peaks_bed}"
        echo -e "\t@ out_prefix ${nucleoatac_out_prefix}"

        sbatch -J ${job_name} ${JOBSCRIPT} ${dataset} \
                            ${nucleoatac_out_prefix} \
                            ${frag_file} \
                            ${peaks_bed} \
                            ${ref_fasta} \
                            ${hits}
    fi

    sleep 2s

done






# # just run the plotting routines for all datasets
# source ../config.sh
# 
# module load python/3.12.1
# eval "$(conda shell.bash hook)"
# conda activate chrombpnet
# 
# finemo_param="counts_v0.23_a0.8_all"
# export R_LIBS_USER=$RENV_SJ
# 
# for i in ${nucleoatac_dir}/*; do
# 
#   j=$(basename -- $i)
#   out_prefix=${nucleoatac_dir%}/$j/$j
#   hits="${hits_reconciled_dir%/}/${j}/${finemo_param}/hits_unique.reconciled.annotated.tsv.gz"
# 
#   # if [[ -f ${out_prefix}.occpeaks.bed.gz && ! -f ${out_prefix}.dyad_binned_dist.pdf ]]; then
# 
#   echo $j
#   
#   python 09-calc_dyad_dist.py ${out_prefix}.occpeaks.bed.gz ${out_prefix} --max_dist 250 --k_nearest 100
#   Rscript 09-plot_nucleoatac.R ${out_prefix}.occpeaks.bed.gz ${hits} ${out_prefix} ${out_prefix}.dyad_binned_distances.tsv
# 
#   # fi
# 
# done

