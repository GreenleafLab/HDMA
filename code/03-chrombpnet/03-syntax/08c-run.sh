#!/bin/bash
#SBATCH --output=../../logs/03-chrombpnet/03/scratch/%x-%j.out
#SBATCH -p akundaje,wjg,sfgf,biochem
#SBATCH --time=06:00:00
#SBATCH -c 2
#SBATCH --mem=16G

source ../config.sh

# motifs -------------------------------------------------------------------
pattern_types=( "pos_patterns" "neg_patterns" )

for pattern_type in ${pattern_types[@]}; do

  echo "#####################################################"
  echo "### @ ${pattern_type}"
  echo "#####################################################"

  # Define the specific motifs and output directory for the current pattern type
  if [[ "${pattern_type}" == "pos_patterns" ]]; then
    motifs_to_process=( "CTCF" "NFY" "YY1" )
  else # This will be "neg_patterns"
    motifs_to_process=( "ZEB" "HIC" "NFY" "YY1" )
  fi
  
  for motif in "${motifs_to_process[@]}"; do
    echo
    echo "--- @ ${motif} for ${pattern_type} ---"
    
    OUT="${base_dir}/03-syntax/08/"
    RESULTS_FILE="${OUT}/jaccard_results.${motif}.${pattern_type}.txt"
    echo $RESULTS_FILE

    # Check if the output file exists
    if [[ -f "${RESULTS_FILE}" ]]; then
        echo "@ done for $RESULTS_FILE"
        continue
    fi

    bash 08c-bedtools_jaccard.sh ${base_dir}/03-syntax/08/${pattern_type}/${motif} ${motif}.${pattern_type}
    
    echo "@ done."
  done
  
done


# regions -------------------------------------------------------------------
echo "#####################################################"
echo "### @ PEAKS"
echo "#####################################################"

region_types=( "Distal" "Promoter" "Exonic" "Intronic")

for region_type in ${region_types[@]}; do

  echo "--- @ ${region_type} ---"
 
  
  OUT="${base_dir}/03-syntax/scratch/02/"
  RESULTS_FILE="${OUT}/jaccard_results.${region_type}.txt"
  echo $RESULTS_FILE

  # Check if the output file exists
  if [[ -f "${RESULTS_FILE}" ]]; then
      echo "@ already done"
      echo ""
      continue
  fi

  bash 02c-bedtools_jaccard.sh ${base_dir}/03-syntax/scratch/02/peaks/${region_type} ${region_type}
  
  echo "@ done."
 
done

