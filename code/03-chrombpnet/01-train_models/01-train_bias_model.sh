#!/bin/bash
#SBATCH --output=../../logs/03-chrombpnet/01/%x-%j.out
#SBATCH -p akundaje
#SBATCH -t 1-0
#SBATCH -c 4
#SBATCH --mem=60G
#SBATCH -G 1

# Purpose: trains an enzymatic bias model for a specific cluster, with a specific
# bias-threshold-factor, using the chrombpnet bias pipeline command, following
# https://github.com/kundajelab/chrombpnet/wiki/Bias-model-training.
# The specific cluster / threshold to use are set in the
# bias_cluster and min_thresh params. Multiple models are trained and the best
# one (which captures Tn5 sequence preferences but not the motifs of real TFs)
# is used for training all bias-factorized models.



# ENVIRONMENT ------------------------------------------------------------------

# load conda environment
eval "$(conda shell.bash hook)"
conda activate chrombpnet

# load modules
module load cuda/11.2.0
module load cudnn/8.1.1.33
module load system
module load libxml2
module load libxslt
module load perl
module load zlib
module load ghostscript
module load cairo
module load pango # Dependency for chrombpnet's make_html()



# PARAMETERS -------------------------------------------------------------------

# source configuration variables
source ../config.sh

# set parameters
min_thresh=0.4
bias_cluster=Heart_c0
bias_fold="fold_0"
data_type=ATAC

ref_fasta="${refs}/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta"
frag_file=${cluster_frags_dir}/fragments/${bias_cluster}__sorted.tsv
peaks_file=${chrombpnet_peaks_dir}/${bias_cluster}__peaks_bpnet.narrowPeak
negatives_file=${negatives_dir}/${bias_cluster}/${bias_fold}/output_negatives.bed
split_file=${split_dir}/${bias_fold}.json
out_dir="${bias_dir%/}/${bias_cluster}_thresh${min_thresh}/"

mkdir -p ${out_dir}



# RUN --------------------------------------------------------------------------
                                   
echo "Running: chrombpnet bias pipeline --genome ${ref_fasta}"
echo "--input-fragment-file ${frag_file}"
echo "--peaks ${peaks_file}"
echo "--nonpeaks ${negatives_file}"
echo "--chr-fold-path ${split_file}"
echo "--bias-threshold-factor ${min_thresh}"
echo "--output-dir ${out_dir}"
echo "--chrom-sizes ${chromsizes}"
echo "--data-type ${data_type}"
                                 
chrombpnet bias pipeline --genome ${ref_fasta} \
                         --input-fragment-file ${frag_file} \
                         --peaks ${peaks_file} \
                         --nonpeaks ${negatives_file} \
                         --chr-fold-path ${split_file} \
                         --bias-threshold-factor ${min_thresh} \
                         --output-dir ${out_dir} \
                         --chrom-sizes ${chromsizes} \
                         --data-type ${data_type}
                                   
