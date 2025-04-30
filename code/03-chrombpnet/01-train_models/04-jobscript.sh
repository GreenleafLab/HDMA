#!/bin/bash
#SBATCH --output=../../logs/03-chrombpnet/01/04/%x-%j.out
#SBATCH -p akundaje,wjg,biochem,sfgf
#SBATCH --time=2-00:00:00
#SBATCH -c 1
#SBATCH --mem=100G
#SBATCH -G 1
#SBATCH --requeue
#SBATCH --open-mode=append

# NOTE: for jobs that timeout with 2 days,
# switch, set a longer time limit and submit to akundaje partition only
# #SBATCH -p akundaje
# #SBATCH --time=3-00:00:00
# otherwise:
# #SBATCH -p akundaje,owners,gpu
# #SBATCH --time=2-00:00:00

ref_fasta=${1}
peaks_file=${2}
model_file=${3}
out_prefix=${4}
chrom_sizes=${5}

function timestamp {
    # Function to get the current time with the new line character
    # removed 
    
    # current time
    date +"%Y-%m-%d_%H-%M-%S" | tr -d '\n'
}

# load conda environment
eval "$(conda shell.bash hook)"
conda activate chrombpnet

module load cuda/11.2
module load cudnn/8.1

echo "--- $(timestamp): Beginning interpretation ---"

chrombpnet contribs_bw --genome ${ref_fasta} \
                                  --regions ${peaks_file} \
                                  --model-h5 ${model_file} \
                                  --output-prefix ${out_prefix} \
                                  --chrom-sizes ${chrom_sizes} \
				  --profile-or-counts counts
                                  # --profile-or-counts # defaults to both

echo "--- $(timestamp): Completed interpretation ---"
