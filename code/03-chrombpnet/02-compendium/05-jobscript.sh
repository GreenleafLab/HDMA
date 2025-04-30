#!/bin/bash
#SBATCH --output=../../logs/03-chrombpnet/02/05/%x-%j.out
#SBATCH -p akundaje,sfgf,gpu,owners
#SBATCH -t 24:00:00
#SBATCH -c 4
#SBATCH --mem=12G
#SBATCH -G 1
#SBATCH --requeue
#SBATCH --open-mode=append

set -eo pipefail

# get args
celltype=${1}
peaks_bed=${2}
modisco_h5=${3}
shaps_h5=${4}
finemo_out=${5}
finemo_out_nocompo=${6}
shap_type=${7}
alpha=${8}
alpha_nocompo=${9}
motifs_nocompo=${10}

echo -e "\t@ ${celltype}"
echo -e "\t@ ${peaks_bed}"
echo -e "\t@ ${modisco_h5}"
echo -e "\t@ ${shaps_h5}"
echo -e "\t@ ${finemo_out}"
echo -e "\t@ ${finemo_out_nocompo}"
echo -e "\t@ ${alpha}"
echo -e "\t@ ${alpha_nocompo}"
echo -e "\t@ ${motifs_nocompo}"

# load conda environment
# NOTE: this call to conda triggers an unbound variable error, hence we 
# use set -eo pipefail instead of -euo pipefail
eval "$(conda shell.bash hook)"
conda activate finemo

# load modules
module load cuda/11.2
module load cudnn/8.1

# print version for debugging
pip freeze | grep finemo


# RUN FINEMO

i=1
outs=( "$finemo_out" "$finemo_out_nocompo" )
for out in ${outs[@]}; do

  echo ${out}

  finemo_npz="${out%/}/intermediate_inputs.npz"

  if [[ -f "${out%/}/hits.bed.gz" ]]; then

    echo "@ found hits."

  else

    echo "[$(date +"%m/%d/%Y (%r)")] starting ${celltype}"
    echo "@ extract-regions-h5"
    finemo extract-regions-chrombpnet-h5 --h5s ${shaps_h5} --out-path ${finemo_npz} --region-width 1000

    echo "[$(date +"%m/%d/%Y (%r)")]"
    echo "@ call-hits"

    # the first time, use all motifs
    if [ $i -eq 1 ]; then
      finemo call-hits -r ${finemo_npz} -m ${modisco_h5} -p ${peaks_bed} --a ${alpha} -o ${out} -b 200
    # the second time, use non-composite motifs only
    elif [ $i -eq 2 ]; then
      finemo call-hits -r ${finemo_npz} -m ${modisco_h5} -p ${peaks_bed} --a ${alpha_nocompo} -o ${out} -b 200 --motifs-include ${motifs_nocompo}
    fi

    # make report
    finemo report -r ${finemo_npz} -H ${out}/hits.tsv -p ${peaks_bed} -m ${modisco_h5} -o ${out} --no-recall --no-seqlets

    echo "[$(date +"%m/%d/%Y (%r)")]"

    # INDEX HITS -----------------------------------------------------------------
    if [[ -f "${out}/hits.bed" ]]; then

    	echo "@ found hits; renaming and indexing"

    	# bgzip and index hits
    	module load biology samtools

    	bgzip -c ${out}/hits.bed > ${out}/hits.bed.gz
    	tabix -p bed ${out}/hits.bed.gz

    fi

  fi

  ((i++))

done

echo "[$(date +"%m/%d/%Y (%r)")] done"
