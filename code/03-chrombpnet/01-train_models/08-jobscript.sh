#!/bin/bash
#SBATCH --output=../../logs/03-chrombpnet/01/08/%x-%j.out
#SBATCH -p akundaje,owners,gpu,wjg
#SBATCH -t 2:00:00
#SBATCH -c 1
#SBATCH --mem=40G
#SBATCH -G 1
#SBATCH --requeue
#SBATCH --open-mode=append

celltype=${1}
peaks_file=${2}
ref_fasta=${3}
chrom_sizes=${4}
out_prefix=${5}
out_key=${6}
model_fold_0=${7}
model_fold_1=${8}
model_fold_2=${9}
model_fold_3=${10}
model_fold_4=${11}

echo "celltype=${1}"
echo "peaks_file=${2}"
echo "ref_fasta=${3}"
echo "chrom_sizes=${4}"
echo "out_prefix=${5}"
echo "out_key=${6}"
echo "model_fold_0=${7}"
echo "model_fold_1=${8}"
echo "model_fold_2=${9}"
echo "model_fold_3=${10}"
echo "model_fold_4=${11}"


# source configuration variables
source ../config.sh

# load conda environment
eval "$(conda shell.bash hook)"
conda activate chrombpnet

module load cuda/11.2
module load cudnn/8.1

echo "[$(date +"%m/%d/%Y (%r)")] starting ${celltype}"

python ./08-predict_and_avg.py --regions ${peaks_file} \
  --genome ${ref_fasta} \
  --chrom-sizes ${chrom_sizes} \
  --output-prefix ${out_prefix} \
  --output-key ${out_key} \
  --output-bed True \
  --batch-size 16 \
  --chrombpnet-model ${model_fold_0} \
  --chrombpnet-model ${model_fold_1} \
  --chrombpnet-model ${model_fold_2} \
  --chrombpnet-model ${model_fold_3} \
  --chrombpnet-model ${model_fold_4}

echo "[$(date +"%m/%d/%Y (%r)")] done!"
