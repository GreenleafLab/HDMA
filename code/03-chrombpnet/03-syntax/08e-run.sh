#!/bin/bash
#SBATCH --output=../../logs/03-chrombpnet/03/scratch/%x-%j.out
#SBATCH -p akundaje,wjg,sfgf,biochem
#SBATCH --time=06:00:00
#SBATCH -c 2
#SBATCH --mem=16G

source ../config.sh

OUT="${base_dir}/03-syntax/08/"

echo "#####################################################"
echo "### @ positive to negative"
echo "#####################################################"

RESULTS_FILE="${OUT}/jaccard_results.pos_to_neg.txt"
echo $RESULTS_FILE

# check if the output file exists
if [[ -f "${RESULTS_FILE}" ]]; then
    echo "@ done for $RESULTS_FILE"
else
    bash 02e-bedtools_jaccard.sh \
      ${base_dir}/03-syntax/08/patterns/pos_patterns/ \
      ${base_dir}/03-syntax/08/patterns/neg_patterns/ \
      pos_to_neg
fi



echo "#####################################################"
echo "### @ NFY direction switch"
echo "#####################################################"

RESULTS_FILE="${OUT}/jaccard_results.NFY_switch.txt"
echo $RESULTS_FILE

# check if the output file exists
if [[ -f "${RESULTS_FILE}" ]]; then
    echo "@ done for $RESULTS_FILE"
else
    bash 02e-bedtools_jaccard.sh \
      ${base_dir}/03-syntax/08/pos_patterns/NFY \
      ${base_dir}/03-syntax/08/neg_patterns/NFY/ \
      NFY_switch
fi



echo "#####################################################"
echo "### @ YY1 direction switch"
echo "#####################################################"

RESULTS_FILE="${OUT}/jaccard_results.YY1_switch.txt"
echo $RESULTS_FILE

# check if the output file exists
if [[ -f "${RESULTS_FILE}" ]]; then
    echo "@ done for $RESULTS_FILE"
else
    bash 02e-bedtools_jaccard.sh \
      ${base_dir}/03-syntax/08/pos_patterns/YY1 \
      ${base_dir}/03-syntax/08/neg_patterns/YY1/ \
      YY1_switch
fi
