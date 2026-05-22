# Author: Selin Jessa
# Purpose: Given a pair of motifs predicted to have cooperativity, run the ISM
# in every cell type to see if it's context specific.
# 
# Run with:
# $ python -u test_coop_across_celltypes.py --motif '446|SRF_TEAD'

# SETUP ------------------------------------------------------------------------

# misc
import sys
import os
import numpy as np
import random
from tqdm import tqdm
import pandas as pd

# io
import h5py
from pyfaidx import Fasta

# ML
import bpnetlite
from bpnetlite import ChromBPNet
from bpnetlite import BPNet
from bpnetlite.bpnite import CountWrapper
from bpnetlite.attribute import deep_lift_shap

import torch
import tangermeme
from tangermeme.io import extract_loci
from tangermeme.io import read_meme
from tangermeme.utils import random_one_hot
from tangermeme.utils import one_hot_encode
from tangermeme.utils import characters
from tangermeme.utils import reverse_complement
from tangermeme.tools.fimo import fimo
from tangermeme.space import space

from scipy.stats import wilcoxon
from scipy.stats import kstest

import pickle
import argparse
import gc
import logging
import datetime

# plotting
import matplotlib.pyplot as plt
import seaborn; seaborn.set_style('white')
import matplotlib.ticker as ticker
import plotnine as pn
import logomaker
from tangermeme.plot import plot_logo

print(torch.__version__)
print(tangermeme.__version__)
print(bpnetlite.__version__)
print(np.__version__)

# os.environ['CUDA_VISIBLE_DEVICES']='0'
os.environ['CUDA_VISIBLE_DEVICES']='MIG-40f43250-998e-586a-ac37-d6520e92590f'

torch.manual_seed(100)

# editable text in PDFs
# https://stackoverflow.com/a/54111532
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42

sys.path.append("..")
import tangermeme_utils as tanutils
from tangermeme_utils.utils import plot_motif
from tangermeme_utils.eval import summarize_across_folds, mean_across_sequences
from tangermeme_utils.wrappers import ChromBPNetWrapper, reshape_folds

from plotnine import (
	ggplot, aes, geom_boxplot, geom_jitter, theme, xlab, ylab, ggtitle, element_blank, element_text, element_rect, scale_fill_manual, coord_cartesian, geom_hline
)





# GLOBAL VARS -------------------------------------------------------------------------

with open("../../DURGA_DIRS.txt", 'r') as f:
	# read the first line
	PROJ_IN = f.readline().strip()

	# read the second line
	PROJ_OUT = f.readline().strip()
	PROJ_OUT = os.path.join(PROJ_OUT, "03-chrombpnet/")


with open("../../AK_PROJ_DIR.txt", 'r') as f:
    KUNDAJE_DIR = f.readline().strip()


INPUTLEN = 2114
OUTPUTLEN = 1000
MIDPOINT = INPUTLEN // 2
SHIFT_REL_INPUT = np.int64((2114 - 1000) / 2)
SHIFT_REL_INPUT

NARROWPEAK_SCHEMA = ["chr", "start", "end", "1", "2", "3", "4", "5", "6", "summit"]

LOGO_ALPHABET = 'ACGT'

# finemo colorscheme
LOGO_COLORS = {"A": '#109648', "C": '#255C99', "G": '#F7B32B', "T": '#D62839'}

# logomaker colorscheme
LOGO_COLORS2= {
		'G': [1, .65, 0],
		'T': [1, 0, 0],
		'C': [0, 0, 1],
		'A': [0, .5, 0]
	}

CHAR_IGNORE = ['QWERYUIOPSDFHJKLZXVBNM']



def test_gpu():

	torch.cuda.device_count()
	dev = torch.device("cuda") if torch.cuda.is_available() else torch.device("cpu")
	dev

	torch.cuda.is_available()

	t1 = torch.randn(1,2).to(dev)
	t1
	t1.device



def parse_args():
	
	parser = argparse.ArgumentParser()
	
	parser.add_argument("--motif", type=str, default=None, help="If present, only run the test for this motif.")
	parser.add_argument("--downstream-only", action="store_true", default=False,
						help="If True, only run the plotting and downstream routines; no marginalizations are performed. Ignores .done file.")
	parser.add_argument("--overwrite", action="store_true", default=False,
						help="If True, overwrite the previous computations. To simply rerun the downstream plotting/summary routines, use --downstream-only.")
	args = parser.parse_args()
	print(args)

	return args





def main(args):

	# set up logging
	logger = logging.getLogger()
	logger.setLevel(logging.WARNING)
	
	# create console handler
	ch1 = logging.StreamHandler()
	ch1.setLevel(logging.WARNING)
	
	# create formatter
	formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
	
	# add formatter to ch
	ch1.setFormatter(formatter)
	
	# add ch to logger
	logger.addHandler(ch1)

	# mlogger = logging.getLogger('matplotlib')
	# mlogger.setLevel(logging.WARNING)
		
 
	# PATHS ------------------------------------------------------------------------
	genome_fa = os.path.join(KUNDAJE_DIR, "refs/hg38/chrombpnet/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta")
	compiled_modisco_h5_path = os.path.join(PROJ_OUT, "02-compendium/modisco_compiled/modisco_compiled_filtered.h5")
	modisco_compiled_unique = os.path.join(PROJ_OUT, "03-syntax/01/motifs_compiled_unique.tsv")
	
	# output paths
	out = os.path.join(PROJ_OUT, "03-syntax/05a/in_silico_marginalization/")

	os.makedirs(out, exist_ok=True)


	# LOAD DATA --------------------------------------------------------------------
	logger.critical("Loading resources...")

	compo_to_test = pd.read_csv("04b-composites_to_test.tsv", sep="\t")
	compo_results = pd.read_csv(os.path.join(PROJ_OUT, "03-syntax/04c/composite_motif_ISM_results.tsv"), sep="\t")

	# check if the motif is in the list
	assert args.motif in compo_to_test.motif_name.values, f"Motif {args.motif} not found in the list of motifs to test."

	comp_motif_name = args.motif

	compo_to_test = compo_to_test[compo_to_test.motif_name == args.motif]
	logger.critical(f"Running with {comp_motif_name}")

	# load list of cell types with models passing QC
	clusters = pd.read_csv(os.path.join(PROJ_OUT, "01-models/qc/chrombpnet_models_keep2.tsv"), sep="\t", header=None)

	# get the first column as a list
	clusters = clusters[0].tolist()
	clusters

	comp_motif_name_safe = comp_motif_name.replace("|", ".").replace("/", ".")
	logger.critical(comp_motif_name_safe)

	# make out directory and paths
	out_dir = os.path.join(out, comp_motif_name_safe)
	os.makedirs(out_dir, exist_ok=True)

	# RUN COOPERATIVITY TESTS ------------------------------------------------------

	# for each composite motif, run in silico marginalziations at optimal spacing
	for cluster in clusters:

		# datetimenow = datetime.datetime.now().strftime("%Y-%m-%d-%H:%M:%S")
		logger.critical(f"--- {cluster}: {comp_motif_name} ---")
		path_done = os.path.join(out_dir, "." + cluster + "__done")

		# setup --------------------------------------------------------------------
		# get info for the experiment

		try:

			# check if done file exists
			if not args.overwrite and not args.downstream_only and os.path.exists(path_done):

				logger.critical(f"{cluster} {comp_motif_name} already done, skipping")
				continue

			else:

				# write a message if done file exists
				if os.path.exists(path_done) and args.downstream_only:
					logger.critical(f"{cluster} {comp_motif_name} is complete; overwriting downstream plots and summaries.")
				elif os.path.exists(path_done):
					logger.critical(f"{cluster} {comp_motif_name} is complete; overwriting computations, plots, and summaries.")

				path_preds_pkl = os.path.join(out_dir, cluster + "__predictions.pkl")
				path_results_df = os.path.join(out_dir, cluster + "__results.tsv")
				path_config_df = os.path.join(out_dir, cluster + "__config.tsv")


				logger.critical(f"{comp_motif_name}")

				# get results for specific composite motif
				row_results = compo_results[compo_results.motif_name == comp_motif_name]
				row_results

				row_test = compo_to_test[compo_to_test.motif_name == comp_motif_name]
				row_test

				seqA = row_results.seqA.values[0]
				seqB = row_results.seqB.values[0]
				comp_type = row_results.category.values[0]
				best_spacing = row_results.best_spacing.values[0]
				best_orientation = row_results.best_orientation.values[0]

				print(best_spacing)
				print(best_orientation)

				motifs = [seqA, seqB]

				logger.critical(motifs)

				# load data and model

				# if the .done file exists, no more computation is needed;
				# don't load the model/backgrounds 
				if not args.downstream_only:

					logger.critical("\tLoading model...")
				
					# load model with five folds

					models = [
						ChromBPNetWrapper(BPNet.from_chrombpnet(filename=os.path.join(chrombpnet_dir, f"01-models/models/bias_Heart_c0_thresh0.4/{cluster}/fold_0/models/chrombpnet_nobias.h5"))),
						ChromBPNetWrapper(BPNet.from_chrombpnet(filename=os.path.join(chrombpnet_dir, f"01-models/models/bias_Heart_c0_thresh0.4/{cluster}/fold_1/models/chrombpnet_nobias.h5"))),
						ChromBPNetWrapper(BPNet.from_chrombpnet(filename=os.path.join(chrombpnet_dir, f"01-models/models/bias_Heart_c0_thresh0.4/{cluster}/fold_2/models/chrombpnet_nobias.h5"))),
						ChromBPNetWrapper(BPNet.from_chrombpnet(filename=os.path.join(chrombpnet_dir, f"01-models/models/bias_Heart_c0_thresh0.4/{cluster}/fold_3/models/chrombpnet_nobias.h5"))),
						ChromBPNetWrapper(BPNet.from_chrombpnet(filename=os.path.join(chrombpnet_dir, f"01-models/models/bias_Heart_c0_thresh0.4/{cluster}/fold_4/models/chrombpnet_nobias.h5")))
					]

					logger.critical("\tLoading backgrounds...")
					# load GC-matched background sequences
					nonpeaks_bed = os.path.join(chrombpnet_dir, f"00-inputs/gc_matched_negatives/{cluster}/fold_0/output_negatives.bed")

					X_nonpeaks = tangermeme.io.extract_loci(nonpeaks_bed, genome_fa, ignore = CHAR_IGNORE)
					logger.critical(X_nonpeaks.shape)
					logger.critical(X_nonpeaks[0])

					# choose 100 random indices of the X_nonpeaks sequences
					rand_indices = np.random.choice(X_nonpeaks.shape[0], 100, replace=False)
					rand_indices

					# ensure tensors are correctly typed and consistent before the function call
					X_nonpeaks = X_nonpeaks.float()
					X_nonpeaks_rand = X_nonpeaks[rand_indices].float()  # convert to float type
					logger.critical(X_nonpeaks_rand.shape)

					# make a config dict
					# write a DF with config info
					config_df = pd.DataFrame({
						"motif_name": [comp_motif_name],
						"category": [comp_type],
						"cluster": [cluster],
						"seqA": [seqA],
						"seqB": [seqB],
						"rand_indices": [",".join(map(str, rand_indices))]
					})

					config_df.transpose()

				else:

					logger.critical("\tRunning downstream plotting & summary routines only...")
					config_df = pd.read_csv(path_config_df, sep="\t")


				# ISM ----------------------------------------------------------------------

				# Spacing -------------------------------------------------------------------
				# check if the spacing predictions have already been computed
				if not args.overwrite and os.path.exists(path_preds_pkl):

					logger.critical("\tLoading precomputed  ISM...")
					with open(path_preds_pkl, "rb") as f:
						preds = pickle.load(f)

				else:

					logger.critical("\tRunning ISM...")

					# the marginalization for the best spacing/orientation
					y_before_list = []
					y_after_list = []

					for model in models:
						y_before, y_after = space(model, X_nonpeaks_rand, motifs, spacing=[[best_spacing]])
						y_before_list.append(y_before)
						y_after_list.append(y_after)

					# stack predictions per fold to get single tensors
					preds = {
						"y_before": reshape_folds(y_before_list),
						"y_after": reshape_folds(y_after_list)
					}
					
					# save the spacing predictions
					with open(path_preds_pkl, "wb") as f:
						pickle.dump(preds, f)

				# calculate the mean counts across sequences, and then mean/sd across folds
				y_after_mean, y_after_sd = tanutils.eval.summarize_across_folds(tanutils.eval.mean_across_sequences(preds['y_after'][1]))

				# calculate the effect in ln(counts) units
				y_effect = preds['y_after'][1] - preds['y_before'][1]

				# calculate mean across sequences, then the mean/sd across folds
				y_effect_mean, y_effect_sd = tanutils.eval.summarize_across_folds(tanutils.eval.mean_across_sequences(y_effect))

				y_after_mean.shape
				y_effect_mean.shape

				y_after_profile_means, y_after_profile_sd = summarize_across_folds(mean_across_sequences(preds['y_after'][0]))

				y_after_profile_means.shape

				# create a dataframe
				results = pd.DataFrame({
					"motif_name": [comp_motif_name],
					"cluster": [cluster],
					"orientation": [best_orientation],
					"spacing": [best_spacing],
					"counts_after_mean": y_after_mean.squeeze().item(),
					"counts_after_sd": y_after_sd.squeeze().item(),
					"effect_mean": y_effect_mean.squeeze().item(),
					"effect_sd": y_effect_sd.squeeze().item()
				})

				results

				# save the outputs ---------------------------------------------------------

				# save the config if it doesn't exist
				if not os.path.exists(path_config_df):
					config_df.to_csv(path_config_df, sep="\t", index=False)
				
				# save the results
				results.to_csv(path_results_df, sep="\t", index=False)

				# make an empty file ".done"
				open(path_done, 'a').close()

				logger.critical("\tDone.")

		except Exception as e:
			logger.critical(f"Error processing {comp_motif_name}: {e}")
			continue

	logger.critical("Done predictions.")
	logger.critical("Aggregating results...")

	preds_after = []
	preds_before = []
	results = []
	clusters_keep = []

	for cluster in clusters:

		path_preds_pkl = os.path.join(out_dir, cluster + "__predictions.pkl")
		path_results_df = os.path.join(out_dir, cluster + "__results.tsv")

		if not os.path.exists(path_preds_pkl):
			print(f"Skipping {cluster}")
			continue

		preds_list = pd.read_pickle(path_preds_pkl)

		y_before_profile_means, _ = summarize_across_folds(mean_across_sequences(preds_list['y_before'][0]))
		y_after_profile_means, _ = summarize_across_folds(mean_across_sequences(preds_list['y_after'][0]))

		preds_before.append(y_before_profile_means.squeeze(0))
		preds_after.append(y_after_profile_means.squeeze(0))

		results_df = pd.read_csv(path_results_df, sep="\t")
		results.append(results_df)
		clusters_keep.append(cluster)

	# stack the predictions
	preds_after = np.stack(preds_after).squeeze(1)
	preds_before = np.stack(preds_before).squeeze(1)

	# predictions to dataframe
	preds_after_df = pd.DataFrame(preds_after, columns=range(1000), index=clusters_keep)
	preds_after_df.reset_index(inplace=True)
	preds_after_df.head()
	preds_after_df.to_csv(os.path.join(out_dir, "aggregated_mean_profile_predictions_after.tsv"), sep="\t", index=False)

	# predictions to dataframe
	preds_before_df = pd.DataFrame(preds_before, columns=range(1000), index=clusters_keep)
	preds_before_df.reset_index(inplace=True)
	preds_before_df.head()
	preds_before_df.to_csv(os.path.join(out_dir, "aggregated_mean_profile_predictions_before.tsv"), sep="\t", index=False)

	# difference
	preds_diff = preds_after - preds_before
	preds_diff_df = pd.DataFrame(preds_diff, columns=range(1000), index=clusters_keep)
	preds_diff_df.reset_index(inplace=True)
	preds_diff_df.head()
	preds_diff_df.to_csv(os.path.join(out_dir, "aggregated_mean_profile_predictions_diff.tsv"), sep="\t", index=False)

	# concatenate the results
	results = pd.concat(results, ignore_index=True)
	results.to_csv(os.path.join(out_dir, "aggregated_results.tsv"), sep="\t", index=False)
				
	logger.critical("Done.")



if __name__ == '__main__':
	args = parse_args()
	main(args)

