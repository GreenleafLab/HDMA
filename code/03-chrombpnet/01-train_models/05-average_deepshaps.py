# Purpose: given the DeepLIFT contribution scores for five folds of a model,
# here we compute the mean across folds per nucleotide and save the output
# as an h5 file.

import deepdish
import numpy as np
import os
import sys

base_dir = sys.argv[1] # expected to contain fold{0,1,2,3,4}
folds    = sys.argv[2] # comma-separated list of folds to keep, e.g. "fold_0,fold_1,fold_2,fold_3,fold_4" 
shaptype = sys.argv[3] # "counts" or "profile"

print(f"averaging {shaptype} deepshaps for [{folds}] in {base_dir}...")

avg_projected_shap = None
raw = None
avg_shap = None

# split folds into a list
folds = folds.split(",")

# loop over folds
for fold in folds:
	print(f"\{fold}")
	fold_dir = f"{base_dir}/{fold}"
	
	# load shap scores
	deepshap = deepdish.io.load(f"{fold_dir}/peaks_shap.{shaptype}_scores.h5")
	
	# sum up the shap scores
	if raw is None:
		avg_projected_shap = deepshap["projected_shap"]["seq"]
		raw = deepshap["raw"]["seq"]
		avg_shap = deepshap["shap"]["seq"]
	else:
		assert(np.array_equal(raw, deepshap["raw"]["seq"]))
		avg_projected_shap += deepshap["projected_shap"]["seq"]
		avg_shap += deepshap["shap"]["seq"]
		
# divide by 5 to average over folds
print("\taveraging...")
avg_projected_shap /= 5
avg_shap /= 5

print("\tsaving...")
deepshap_output = {"projected_shap": {"seq": avg_projected_shap},
                   "raw": {"seq": raw},
                   "shap": {"seq": avg_shap}}
deepdish.io.save(f"{base_dir}/average_shaps.{shaptype}.h5", deepshap_output)

# done.
