#!/usr/bin/bash
#SBATCH --output=../../logs/03-chrombpnet/02/03/%x-%j.out
#SBATCH -p akundaje,wjg,biochem,sfgf
#SBATCH -t 12:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=80G

set -euo pipefail

# N.B. regarding resources: merging < 80 or so patterns only takes 2-3 hours.
# Merging 1000 patterns takes up to 2 days.

output_dir=${1}
model_head=${2}
cluster_key=${3}
modisco_dir=${4}
contribs_dir=${5}
batch=${6}

# load conda env
eval "$(conda shell.bash hook)"
conda activate modiscolite

module load system
module load cairo

echo "-----------------"
echo "starting"

echo "@ got args: "
echo "--out-dir ${output_dir}"
echo "--model-head ${model_head}"
echo "--cluster-key ${cluster_key}"
echo "--modisco-dir ${modisco_dir}"
echo "--contribs-dir ${contribs_dir}"
echo "--batch ${batch}"

# use -u unbuffered option to get live printed output to follow progress via logs
# https://stackoverflow.com/questions/33178514/how-do-i-save-print-statements-when-running-a-program-in-slurm
python -u 03-merge_modisco.py --out-dir ${output_dir} \
                    --model-head ${model_head} \
                    --cluster-key ${cluster_key} \
                    --modisco-dir ${modisco_dir} \
                    --contribs-dir ${contribs_dir} \
                    --batch ${batch}
                    
echo "done!"
echo "-----------------"