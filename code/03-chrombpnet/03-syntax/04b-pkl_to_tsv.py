# Helper script to extract profile predictions for a given composite motif
# at some specific orientation, for the first 10 spacings, summarized across all 100 sequences
# and 5 folds, and save to a TSV for plotting in R.

import sys
import os
import pickle
import numpy as np
import pandas as pd
import argparse

sys.path.append("..")
import tangermeme_utils as tanutils
from tangermeme_utils.eval import summarize_across_folds, mean_across_sequences

# set up argparsing
parser = argparse.ArgumentParser(description='Extract predicted profiles to TSV')
parser.add_argument('--motif', type=str, help='Motif name')
parser.add_argument('--orientation', type=str, help='Orientation to extract')
args = parser.parse_args()

motif_safe = args.motif.replace("/", ".")

with open("../../ROOT_DIR.txt", 'r') as f:
    PROJ_DIR = f.readline().strip()

with open("../../AK_PROJ_DIR.txt", 'r') as f:
    KUNDAJE_DIR = f.readline().strip()

out_path = os.path.join(PROJ_DIR, "output/03-chrombpnet/03-syntax/04b/in_silico_marginalization/", motif_safe)
print(out_path)

# load the predictions across arrangements
pkl_path = os.path.join(out_path, "predictions.pkl")

with open(pkl_path, 'rb') as f:
	spacing_preds = pickle.load(f)

# extract the orientation
subsetted_preds = spacing_preds[args.orientation]['y_after'][0][:, :, 0:10, :, :]
print(subsetted_preds.shape)
print(type(subsetted_preds))
profile_mean, profile_sd = summarize_across_folds(mean_across_sequences(subsetted_preds))

spacing_list = []
position_list = []
mean_list = []
sd_list = []

for spacing_idx in range(profile_mean.shape[0]):
	for position_idx in range(profile_mean.shape[2]):
		spacing_list.append(spacing_idx)
		position_list.append(position_idx)
		mean_list.append(profile_mean[spacing_idx, 0, position_idx].item())
		sd_list.append(profile_sd[spacing_idx, 0, position_idx].item())

# make a df
data = {
	"spacing": spacing_list,
	"position": position_list,
	"mean": mean_list,
	"sd": sd_list,
}

df = pd.DataFrame(data)
df.to_csv(os.path.join(out_path, f"profile_preds_{args.orientation}.tsv"), sep="\t", index=False)
