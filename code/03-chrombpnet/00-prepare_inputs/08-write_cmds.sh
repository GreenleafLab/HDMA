#!/usr/bin/bash

# Purpose: generate nonpeaks (GC-matched negatives), which will be used as background
# regions during model training. This script writes the commands, one per fold
# per cell type.
# Info: https://github.com/kundajelab/chrombpnet/wiki/Preprocessing#generate-non-peaks-background-regions

# SETUP ------------------------------------------------------------------------

# load conda environment
eval "$(conda shell.bash hook)"
conda activate chrombpnet

# source configuration variables
source ../config.sh

ref_fasta="${refs}/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta"

input_len=2114
stride=1000


# CONSTRUCT COMMANDS ------------------------------------------------

# get all clusters
datasets=$(ls ${peaks_dir%/}/*__peaks_overlap_filtered.narrowPeak | xargs -n 1 -I {} basename {} __peaks_overlap_filtered.narrowPeak)	

cmdfile=08-commands.sh
rm $cmdfile
touch $cmdfile

for cluster in ${datasets[@]}; do
  
  echo ${cluster}
    
  peak_file="${chrombpnet_peaks_dir}/${cluster}__peaks_bpnet.narrowPeak"
  
  cluster_dir=${negatives_dir%/}/${cluster}
  [[ -d ${cluster_dir} ]] || mkdir -p ${cluster_dir}
  
  for fold in {0..4}; do
    fold_name=fold_${fold}
    echo -e "\t${fold_name}"
  
    split_file=${split_dir%/}/${fold_name}.json
  
    fold_dir=${cluster_dir%/}/${fold_name}
    [[ -d ${fold_dir} ]] || mkdir -p ${fold_dir}
    fold_out="${fold_dir}/output"
  
    job_name="08-neg_${cluster}_${fold_name}"
    JOBSCRIPT=08-get_negatives.sh
  
    # submit jobs
    done_file="${fold_dir}/output_negatives.bed"
    
    if [[ ! -f "$done_file" ]]; then
      rm -rf "${fold_dir}/output_auxiliary"
      echo "${JOBSCRIPT} ${ref_fasta} ${chromsizes} ${blacklist} ${peak_file} ${input_len} ${stride} ${fold_out} ${split_file}" >> $cmdfile
    else
      echo "${dataset} ${fold} done; skipping"
    fi
    
  done

done
