#!/bin/bash

parent_dir=../../output/04-enhancers/03/results/
mkdir -p $parent_dir/staging

find $parent_dir -type f -path '*/Predictions/EnhancerPredictionsFull_threshold0.013_self_promoter.tsv' | while read -r file; do
  subdir=$(echo "$file" | awk -F'/' '{print $(NF-2)}')
  echo $subdir
  cp "$file" "$parent_dir/staging/${subdir}_EnhancerPredictionsFull_threshold0.013_self_promoter.tsv"
done

cd $parent_dir/staging
zip ../abc_enhancer_predictions_threshold0.013.zip *
