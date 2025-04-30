# Author: Selin Jessa
# Purpose: Test cooperativity between motifs in a composite motif using in silico marginalization.
# This script loops over pairs of motifs and their high-affinity/representative sequences,
# and runs the spacing experiment for the two motifs across all possible orientations
# and 0-200 bp of distance between motifs. These results are compared to the ISM
# for each motif alone. The sequences with the optimal syntax for the two motifs
# are interpreted (basically giving a marginal DeepSHAP) to confirm that the inserted
# motifs drive the predicted accessibility. Finally this script has several
# output and diagnostic plots for the results.
#
# This code relies heavily on the tangermeme package by Jacob Schreiber.
#
# Run on durga on Kundaje lab cluster, with the following command:
# $ python -u test_cooperativity.py
#
# Run only plotting / downstream routines for all motifs with:
# $ python -u test_cooperativity.py --downstream-only
#
# Debug for one motif with:
# $ python -u test_cooperativity.py --motif '434|SOX_SOX#1'



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
from bpnetlite.bpnet import CountWrapper
from bpnetlite.attribute import deep_lift_shap

import torch
import tangermeme
from tangermeme.io import extract_loci
from tangermeme.io import read_meme
from tangermeme.utils import random_one_hot
from tangermeme.utils import one_hot_encode
from tangermeme.utils import characters
from tangermeme.utils import reverse_complement
from tangermeme.utils import _validate_input
from tangermeme.utils import _cast_as_tensor
from tangermeme.tools.fimo import fimo
from tangermeme.space import space
from tangermeme.ersatz import multisubstitute
from tangermeme.ersatz import substitute
from tangermeme.marginalize import marginalize

from scipy.stats import wilcoxon
from scipy.stats import kstest
from scipy.stats import entropy

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





# CONSTANTS -------------------------------------------------------------------------

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


# TESTING ------------------------------------------------------------------------

def test_gpu():

	torch.cuda.device_count()
	dev = torch.device("cuda") if torch.cuda.is_available() else torch.device("cpu")
	dev

	torch.cuda.is_available()

	t1 = torch.randn(1,2).to(dev)
	t1
	t1.device



# ARGUMENTS ------------------------------------------------------------------------

def parse_args():
	
	parser = argparse.ArgumentParser()
	
	parser.add_argument("--motif", type=str, default=None, help="If present, only run the test for this motif.")
	parser.add_argument("--downstream-only", action="store_true", default=False,
						help="If True, only run the plotting and downstream routines; no marginalizations are performed. Ignores .done file.")
	parser.add_argument("--overwrite", action="store_true", default=False,
						help="If True, overwrite the previous computations. To simply rerun the downstream plotting/summary routines, use --downstream-only.")
	parser.add_argument("--overwrite-interpretation", action="store_true", default=False,
						help="If True, overwrite the previous interpretation. To simply rerun the downstream plotting/summary routines, use --downstream-only.")
	parser.add_argument("--interpret", action="store_true", default=False,
						help="If True, run the interpretation routines. Default: False.")
	parser.add_argument("--interpret-only", action="store_true", default=False,
						help="If True, only run the interpretation routines (and downstream plotting). Default: False.")
	parser.add_argument("--n-spaces", type=int, default=201, help="The number of spacings to test. Default: 201, which will test spacings from 0 to 200 bp.")
	args = parser.parse_args()
	print(args)

	return args



# SPACING EXPERIMENT ----------------------------------------------------------------

def test_spacing(models, motifs, seqs, key, n_spaces = 20):
	"""
	A helper function to run the spacing experiment for one motif pair
	and save the outputs as dataframe.

	Parameters
	----------
	models: List of models as output by ChromBPNetWrapper
		The model to use for the experiment.
	motifs: list of str, expecting two motifs
		The motifs to test for cooperativity.
	seqs: torch.Tensor
		The background sequences to test the motifs on.
	key: str
		The key to use for the dataframe to label the experiment.
	n_spaces: int
		The number of spacings to test, default is 20.

	Returns
	-------
	y_before: torch.Tensor
		The predictions before marginalization.
	y_after: torch.Tensor
		The predictions after marginalization.
	df: pd.DataFrame
		The summarized effects per spacing.
	"""

	# the marginalization per spacing
	y_before_list = []
	y_after_list = []
	for model in models:
		y_before, y_after = space(model, seqs, motifs, spacing=torch.arange(n_spaces)[:, None])
		y_before_list.append(y_before)
		y_after_list.append(y_after)

	# stack predictions by output
	y_before = reshape_folds(y_before_list)
	y_after = reshape_folds(y_after_list)

	# calculate the mean counts across sequences, and then mean/sd across folds
	y_after_mean, y_after_sd = tanutils.eval.summarize_across_folds(tanutils.eval.mean_across_sequences(y_after[1]))

	# calculate the effect in ln(counts) units
	y_effect = y_after[1] - y_before[1]

	# calculate mean across sequences, then the mean/sd across folds
	y_effect_mean, y_effect_sd = tanutils.eval.summarize_across_folds(tanutils.eval.mean_across_sequences(y_effect))

	# create a dataframe
	df = pd.DataFrame({
		"key": [key] * n_spaces,
		"spacing": np.arange(n_spaces),
		"counts_after_mean": y_after_mean.squeeze(),
		"counts_after_sd": y_after_sd.squeeze(),
		"effect_mean": y_effect_mean.squeeze(),
		"effect_sd": y_effect_sd.squeeze()
	})

	return y_before, y_after, df




# CALCULATIONS ---------------------------------------------------------------------

def is_palindrome(seq):
	"""
	Check if a sequence is a palindrome.
	"""

	return seq == reverse_complement(seq)



def calc_effect_Z_score(effects):
	"""
	Calculate the Z-scores for the marginal effects of the motif pairs across
	spacings and orientations.

	Z = (X - mean) / std
	"""

	mean = np.mean(effects)
	std = np.std(effects)
	z = (effects - mean) / std

	return z



def calc_effect_entropy(effects):
	"""
	Calculate the entropy of the effects for the motif pairs across
	spacings and orientations.

	We expect *lower* entropy when the uncertainty is low, meaning there is a
	specific orientation and spacing that is favored.

	We expect a *higher* entropy when the uncertainty is high, because the
	distribution is more uniform.
	"""

	# convert log fold changes to probabilities so that they sum to 1

	# shift in case any are below 0
	effects_shifted = effects - np.min(effects)

	# normalize to convert to probabilities
	probabilities = effects_shifted / np.sum(effects_shifted)
		
	# compute entropy
	S = entropy(probabilities)

	return S





# PLOTTING ------------------------------------------------------------------------

def plot_effect_spacing(spacing_dfs, input_motifs, fig_path):
	"""
	Plot the marginal effects for the motif pairs at every spacing in every 
	orientation during in silico marginalization, as a bar plot.
	"""

	n_col = len(input_motifs)
	fig, axes = plt.subplots(n_col, 1, figsize=(10, n_col * 3))
	plt.set_loglevel('WARNING') 

	# iterate over the motif pairs and plot
	for i, key in enumerate(input_motifs.keys()):

		# get the effects from the dataframe
		effects = spacing_dfs[spacing_dfs.key == key].effect_mean.values
		effects_sd = spacing_dfs[spacing_dfs.key == key].effect_sd.values
		
		# plot bar and errorbars
		axes[i].bar(range(len(effects)), effects, alpha=1, label=key, color="red")
		axes[i].errorbar(range(len(effects)), effects, yerr=effects_sd, fmt='none', ecolor='black', elinewidth=1, capsize=0)
		
		# set titles and labels
		axes[i].set_title(f"{key}: {input_motifs[key][0]}, {input_motifs[key][1]}")
		axes[i].set_xlabel("Spacing between motifs (bp)", fontsize=10)
		axes[i].xaxis.set_major_locator(ticker.MaxNLocator(integer=True, nbins=5))
		axes[i].xaxis.set_minor_locator(ticker.AutoMinorLocator())
		axes[i].tick_params(axis='x', which='both', length=2, direction='out', top=False, bottom=True)
		axes[i].set_xlim([0, args.n_spaces])
		axes[i].set_ylim([0, 3])
		# axes[i].legend()

	axes[0].set_ylabel("Mean marginal effect (ln counts)")

	# joint title
	plt.suptitle("Marginal effects for motif pairs during in silico marginalization", fontsize=12)
	plt.tight_layout()
	plt.savefig(fig_path, dpi=300);



def plot_effect_boxplot(df_pred_effect, best_orientation, best_spacing, fig_path):
	"""
	Plot the marginale effects for the motif pairs at the best spacing and orientation,
	as well as the fold changes for the solo effects and sum of independent effects.
	"""

	df_pred_effect_long = df_pred_effect.melt(id_vars = "example")
	df_pred_effect_long.head()	

	# boxplot of predicted log counts
	# set order of the variables
	df_pred_effect_long["variable"] = pd.Categorical(
		values = df_pred_effect_long["variable"],
		categories = ['Joint', 'A + B', 'A', 'B'],
		ordered = True)

	p = (	
		ggplot(df_pred_effect_long, aes(x = "variable", y = "value"))
		+ geom_boxplot(aes(fill = "variable"))
		+ ggtitle(f"In silico marginalizations for motifs individually and togther \nBest spacing/orientation: {best_orientation} @ {best_spacing} bp")
		+ ylab("Predicted marginal effect (ln counts)")
		+ xlab("Motifs")
		+ scale_fill_manual(values = {"Joint": 'red',
									"A + B": '#741fb5',
									"A": "blue",
									"B": '#4e54ed'})
		+ geom_hline(yintercept=0, linetype="dashed", color = "black")
		+ theme(axis_text_x=element_text(angle = 45),
				panel_background=element_rect(fill="white"),
				panel_border=element_rect(color = "black", size = 0.5),
				figure_size=(6, 4))
		)

	p.save(filename = fig_path, width=6, height=4, units='in', dpi=300, verbose=False);



def plot_effect_zscores(df_z, input_motifs, fig_path):
	"""
	Plot Z-scored marginal effects for motif pairs across spacings and orientations.
	"""

	palette_blue = seaborn.color_palette("Blues", len(input_motifs)).as_hex()
	palette = { key: palette_blue[i] for i, key in enumerate(input_motifs.keys()) }

	# set order of the levels in orientation as the keys of the input motifs
	df_z["orientation"] = pd.Categorical(df_z["orientation"], categories = list(input_motifs.keys()), ordered = True)

	p = ( ggplot(df_z, aes(x = "orientation", y = "Z"))
		+ geom_jitter(aes(fill = "orientation"), alpha = 0.6, size = 4, width = 0.1)
		+ ggtitle("Z-scores for marginal effects across \nspacings/orientations during ISM")
		+ ylab("Z-score")
		+ xlab("Orientation")
		+ scale_fill_manual(values = palette)
		+ geom_hline(yintercept=[-3, -2, 2, 3], linetype="dashed", color = "black")
		+ theme(axis_text_x=element_text(angle = 45),
				panel_background=element_rect(fill="white"),
				panel_border=element_rect(color = "black", size = 0.5),
				figure_size=(6, 4))
	)

	p.save(filename = fig_path, width=6, height=4, units='in', dpi=300, verbose=False);



def plot_profiles_spacing(profile_preds_summarized, input_motifs, fig_path):
	"""
	Plot the predicted profile for the background and edited sequence for insertions of
	motif pairs at every orientation and the first 10 bp of spacings.
	"""

	range_idx = range(250, 750)

	# get max of profiles across all orientations
	maxes = []
	for key in input_motifs.keys():

		max_per_key = torch.max(profile_preds_summarized[key]['y_after_means'])
		maxes.append(max_per_key)

	y_max = np.max([np.max(maxes), 1])

	# make ten small subplots in a row, for as many rows as input motifs
	fig, ax = plt.subplots(len(input_motifs), 10, figsize=(18, len(input_motifs)*2.5))
	plt.set_loglevel('WARNING') 

	# i indexes the motif pair
	for i, key in enumerate(input_motifs.keys()):
		
		y_before_mean = profile_preds_summarized[key]['y_before_means']
		y_before_sd = profile_preds_summarized[key]['y_before_sd']
		y_after_mean = profile_preds_summarized[key]['y_after_means']
		y_after_sd = profile_preds_summarized[key]['y_after_sd']

		# get max profiles across all spacings
		# max_profile = torch.max(y_after_mean)
		# y_max = np.max([max_profile, 1])

		# j indexes the spacing
		for j in range(10):

			ax[i, j].plot(y_before_mean[j].T[range_idx], label="background", color="black", alpha=0.2, lw=0.8)
			ax[i, j].fill_between(range(500),
								(y_before_mean[j].T[range_idx] - y_before_sd[j].T[range_idx]).squeeze(),
								(y_before_mean[j].T[range_idx] + y_before_sd[j].T[range_idx]).squeeze(),
								color="black", alpha=0.1)

			ax[i, j].plot(y_after_mean[j].T[range_idx], label=f"edited \n{key}", color="red", alpha=0.8, lw=0.8)
			ax[i, j].fill_between(range(500),
								(y_after_mean[j].T[range_idx] - y_after_sd[j].T[range_idx]).squeeze(),
								(y_after_mean[j].T[range_idx] + y_after_sd[j].T[range_idx]).squeeze(),
								color="red", alpha=0.1)
			ax[i, j].set_ylim([0, y_max])

			# choose labels for the x-axis, only if it's the last row of plots
			if i!=(len(input_motifs)-1):
				ax[i, j].set_xticks([])
			elif i==(len(input_motifs)-1):
				ax[i, j].set_xticks(ticks = [0, 250, 500], labels = [-250, 0, 250])

			# column titles, if it's the first row of plots
			if i == 0:
				ax[i, j].set_title(f"{j} bp")

			# remove labels for the y-axis, if it's not the first column of plots
			if j != 0:
				ax[i, j].set_yticks([])
				ax[i, j].set_yticks([])

		ax[i, 0].set_ylabel("Avg. predicted profile")
		ax[i, 9].legend()

	# shared title
	fig.suptitle("Predicted profiles for motifs inserted at different spacings", fontsize=12)

	# shared x-axis label
	fig.text(0.5, 0.01, "Position relative to insertion (bp)", ha='center', fontsize=10)

	plt.savefig(fig_path, dpi=300);



def plot_profiles_solo(solo_profile_preds_summarized, profile_preds_summarized, best_orientation, best_spacing, fig_path):
	"""
	Plot the predicted profiles for the best spacing and orientation, as well as the independent effects
	and sum of independent effects.
	"""

	plt.figure(figsize=(7, 4))
	plt.set_loglevel('WARNING') 

	# plot the background from one experiment
	plt.plot(solo_profile_preds_summarized['A']['y_before_means'].T, label="background", color="lightgray", alpha=0.8, lw=0.8)
	plt.fill_between(range(1000),
	(solo_profile_preds_summarized['A']['y_before_means'].T - solo_profile_preds_summarized['A']['y_before_sd'].T).squeeze(),
	(solo_profile_preds_summarized['A']['y_before_means'].T + solo_profile_preds_summarized['A']['y_before_sd'].T).squeeze(), color="lightgray", alpha=0.1)

	# get max profiles across all spacings
	max_profile = torch.max(profile_preds_summarized[best_orientation]['y_after_means'][best_spacing])
	y_max = np.max([max_profile, 1])


	# plot the solo effects
	for key, color in zip(['A', 'B', 'A + B'], ['#4e54ed', 'blue', '#741fb5']):
		plt.plot(solo_profile_preds_summarized[key]['y_after_means'].T, label=key, color=color, alpha=0.8, lw=0.8)
		plt.fill_between(range(1000), 
						(solo_profile_preds_summarized[key]['y_after_means'].T - solo_profile_preds_summarized[key]['y_after_sd'].T).squeeze(), 
						(solo_profile_preds_summarized[key]['y_after_means'].T + solo_profile_preds_summarized[key]['y_after_sd'].T).squeeze(), 
						color=color, alpha=0.1)

	# plot the profile for the best spacing/orientation
	plt.plot(profile_preds_summarized[best_orientation]['y_after_means'][best_spacing].T, label=f"{best_orientation} @ {best_spacing} bp", color="red", alpha=0.8, lw=0.8)
	plt.fill_between(range(1000), 
					(profile_preds_summarized[best_orientation]['y_after_means'][best_spacing].T - profile_preds_summarized[best_orientation]['y_after_sd'][best_spacing].T).squeeze(), 
					(profile_preds_summarized[best_orientation]['y_after_means'][best_spacing].T + profile_preds_summarized[best_orientation]['y_after_sd'][best_spacing].T).squeeze(), 
					color="#ff7077", alpha=0.1)

	ax = plt.gca()
	ax.xaxis.set_major_locator(ticker.MaxNLocator(integer=True, nbins = 10))
	ax.xaxis.set_minor_locator(ticker.AutoMinorLocator())
	ax.tick_params(axis='x', which='both', length=5, direction='out', top=False, bottom=True)

	# put legend outside
	plt.legend(loc='upper left', bbox_to_anchor=(1, 1))
	plt.title("Predicted profiles for joint and independent motifs")
	plt.ylim([0, y_max])
	plt.ylabel("Avg. predicted profile")
	plt.xlabel("Position relative to insertion (bp)")
	plt.xticks(ticks = [0, 500, 1000], labels = [-500, 0, 500])
	plt.savefig(fig_path, dpi=300, bbox_inches='tight');



def plot_effect_solo(solo_dfs, best_orientation, best_spacing, fig_path):
	"""
	Plot the effects for the solo effects and the sum of independent effects
	as bar plots.
	"""

	fig, ax = plt.subplots(1, 1, figsize=(5, 3))
	plt.set_loglevel('WARNING') 

	ax.bar(solo_dfs.key, solo_dfs.effect_mean, alpha=1, label=best_orientation, color = ['#4e54ed', 'blue', '#741fb5'])
	ax.errorbar(solo_dfs.key, solo_dfs.effect_mean, yerr=solo_dfs.effect_sd, fmt='none', ecolor='black', elinewidth=1, capsize=0)

	ax.set_ylabel("Mean marginal effects (ln counts)")
	ax.set_ylim([0, 3])
	plt.title("Marginal effects in silico marginalizations of solo motifs", fontsize=12)

	plt.savefig(fig_path, dpi=300);



def plot_effect_spacing_lines(spacing_df, solo_df, motif_name, fig_path, ymax = 1):

	# plot this as a line plot
	fig, axes = plt.subplots(4, 1, figsize=(13, 11))

	n_spaces = 201

	# iterate over the motif pairs and plot
	for i, key in enumerate(spacing_df.key.unique()):

		# get the l2fcs from the dataframe
		effects = spacing_df[spacing_df.key == key].effect_mean.values
		effects_sd = spacing_df[spacing_df.key == key].effect_sd.values
		
		# plot bar and errorbars
		axes[i].plot(range(len(effects)), effects, alpha=1, label=key, color="red")
		axes[i].fill_between(range(len(effects)), effects + effects_sd, effects - effects_sd, color = "red", alpha = 0.1)

		# set titles and labels
		axes[i].set_title(f"{key}")
		axes[i].xaxis.set_major_locator(ticker.MaxNLocator(integer=True, nbins=5))
		axes[i].xaxis.set_minor_locator(ticker.AutoMinorLocator())
		axes[i].tick_params(axis='x', which='both', length=2, direction='out', top=False, bottom=True)
		axes[i].set_xlim([0, n_spaces])
		axes[i].set_ylim([0, ymax])
		# axes[i].legend()

	color_pal = {"Joint": 'red',
	"A + B": '#741fb5',
	"A": "blue",
	"B": '#4e54ed'}

	# plot a straight line
	def plot_effect_with_error(ax, x, y, yerr, label, color, alpha=0.1):
		ax.plot(x, y, alpha=1, label=label, color=color)
		ax.fill_between(x, y + yerr, y - yerr, color=color, alpha=alpha)

	for i in range(4):
		plot_effect_with_error(axes[i], range(n_spaces), np.repeat(solo_df[solo_df.key == "A"].effect_mean.values, n_spaces), np.repeat(solo_df[solo_df.key == "A"].effect_sd.values, n_spaces), "A", color_pal['A'])
		plot_effect_with_error(axes[i], range(n_spaces), np.repeat(solo_df[solo_df.key == "B"].effect_mean.values, n_spaces), np.repeat(solo_df[solo_df.key == "B"].effect_sd.values, n_spaces), "B", color_pal['B'])
		plot_effect_with_error(axes[i], range(n_spaces), np.repeat(solo_df[solo_df.key == "A + B"].effect_mean.values, n_spaces), np.repeat(solo_df[solo_df.key == "A + B"].effect_sd.values, n_spaces), "A + B", color_pal['A + B'])

	axes[1].set_ylabel("Mean marginal effect")
	axes[0].legend()
	axes[3].set_xlabel("Spacing between motifs (bp)", fontsize=10)

	# joint title
	plt.suptitle(f"Mean effects for motif pairs during ISM: {motif_name}", fontsize=12)
	plt.tight_layout()
	plt.savefig(fig_path, dpi = 300);



def plot_mean_contributions(profile_preds_summarized, contributions_mean, coords, input_motifs, best_spacing, best_orientation, fig_path):
	"""
	Plot the mean attributions for the best spacing and orientation. Specifically,
	we plot the average predicted profile across the 1 kb, the average contribution
	scores across the central 1 kb of the input (predicted window), and the zoomed
	in contributions at the center of the predicted window, as a sequence logo.

	This function has some careful handling of starts and end positions, because
	the input window (2,114) and the prediction window (1,000) are not the same.
	So we often want to subset regions _relative_ to one of these windows, for example
	the specific positions where the motifs were inserted.
	"""

	fig, ax = plt.subplots(3, 1, figsize = (11, 6))
	plt.set_loglevel('WARNING') 

	# plto the background
	ax[0].plot(profile_preds_summarized[best_orientation]['y_before_means'][best_spacing].T, linewidth=1, label = "mean prediction (background)", color = 'lightgray')
	ax[0].fill_between(range(1000),
					(profile_preds_summarized[best_orientation]['y_before_means'][best_spacing].T - profile_preds_summarized[best_orientation]['y_before_sd'][best_spacing].T).squeeze(),
					(profile_preds_summarized[best_orientation]['y_before_means'][best_spacing].T + profile_preds_summarized[best_orientation]['y_before_sd'][best_spacing].T).squeeze(),
			color="lightgray", alpha=0.1, label = "sd across folds")

	# plot the predicted profile for the top hit
	ax[0].plot(profile_preds_summarized[best_orientation]['y_after_means'][best_spacing].T, linewidth=1, label = "mean prediction (edited)", color = 'red')
	ax[0].fill_between(range(1000),
					(profile_preds_summarized[best_orientation]['y_after_means'][best_spacing].T - profile_preds_summarized[best_orientation]['y_after_sd'][best_spacing].T).squeeze(),
					(profile_preds_summarized[best_orientation]['y_after_means'][best_spacing].T + profile_preds_summarized[best_orientation]['y_after_sd'][best_spacing].T).squeeze(),
			color="red", alpha=0.1, label = "sd across folds")
	
	# labels and axes
	ax[0].set_ylabel("Predicted profile")
	ax[0].legend(fontsize = 8)
	ax[0].set_xticks(ticks = [0, 100, 200, 300, 400, 500, 600, 700, 800, 900, 999],
					 labels = [-500, -400, -300, -200, -100, 0, +100, +200, +300, +400, +500])

	# semi-transparent rectangle to indicate the zoom at the bottom
	ax[0].axvspan(coords['start_A'] - SHIFT_REL_INPUT, coords['end_A'] - SHIFT_REL_INPUT, alpha=0.4, color='yellow')
	ax[0].axvspan(coords['start_B'] - SHIFT_REL_INPUT, coords['end_B'] - SHIFT_REL_INPUT, alpha=0.4, color='yellow')

	# contribs
	# plot the sum of contributions across nucleotides, per position
	ax[1].plot(contributions_mean.sum(0)[ MIDPOINT - 500 : MIDPOINT + 500 ], linewidth=1, label = "contrib scores", color = "blue")
	ax[1].set_ylabel("Contrib.")
	ax[1].legend(fontsize = 8)
	ax[1].axvspan(coords['start_A'] - SHIFT_REL_INPUT, coords['end_A'] - SHIFT_REL_INPUT, alpha=0.4, color='yellow')
	ax[1].axvspan(coords['start_B'] - SHIFT_REL_INPUT, coords['end_B'] - SHIFT_REL_INPUT, alpha=0.4, color='yellow')
	ax[1].set_xticks(ticks = [0, 100, 200, 300, 400, 500, 600, 700, 800, 900, 999],
				  	 labels = [-500, -400, -300, -200, -100, 0, +100, +200, +300, +400, +500])

	tangermeme.plot.plot_logo(np.float64(contributions_mean), start = coords['zoom_start'], end = coords['zoom_end'], ax=ax[2])
	# adjust coords by 0.5 to include whole nucleotides in the highlight
	ax[2].axvspan(coords['start_A'] - coords['zoom_start'], coords['end_A'] - coords['zoom_start'], alpha=0.3, color='yellow')
	ax[2].axvspan(coords['start_B'] - coords['zoom_start'], coords['end_B'] - coords['zoom_start'], alpha=0.3, color='yellow')
	ax[2].set_ylabel("Contrib.")
	ax[2].set_xticks(
		ticks = [0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 99],
		labels = [-50, -40, -30, -20, -10, 0, 10, 20, 30, 40, 50])

	# put a title at the very top
	plt.suptitle(f"Attributions for {best_orientation} @ {best_spacing} bp: {input_motifs[best_orientation][0]}, {input_motifs[best_orientation][1]}", fontsize = 12)
	plt.xlabel("Genomic position")
	plt.tight_layout()

	plt.savefig(fig_path, bbox_inches='tight', dpi = 100);



def plot_contributions_heatmap(contribs_sum_mean_across_folds, coords, best_orientation, best_spacing, fig_path):
	"""
	Plot a heatmap of the contribution scores in the central 100bp, for each sequence.
	"""

	# make a heatmap
	fig, ax = plt.subplots(1, 1, figsize = (5, 7))
	plt.set_loglevel('WARNING') 
	seaborn.heatmap(contribs_sum_mean_across_folds[:, coords['zoom_start']:coords['zoom_end']], cmap = "magma", ax = ax)

	# highlight regions w/ motifs inserted with lines
	ax.axvline(x = coords['start_A'] - coords['zoom_start'], color = "black", lw = 1, ls = "--")
	ax.axvline(x = coords['end_A'] - coords['zoom_start'], color = "black", lw = 1, ls = "--")
	ax.axvline(x = coords['start_B'] - coords['zoom_start'], color = "black", lw = 1, ls = "--")
	ax.axvline(x = coords['end_B'] - coords['zoom_start'], color = "black", lw = 1, ls = "--")

	plt.ylabel("Sequence")
	plt.xlabel("Genomic position (central 100 bp)")
	plt.yticks([])

	plt.xticks(
		ticks = [0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 99],
		labels = [-50, -40, -30, -20, -10, 0, 10, 20, 30, 40, 50])

	plt.title(f"Contributions for {best_orientation} @ {best_spacing} bp \n(mean across 5 folds per sequence)")
	plt.tight_layout()
	plt.savefig(fig_path, dpi = 300);



def plot_profile_heatmap_z(preds, coords, best_orientation, best_spacing, fig_path):
	"""
	Plot a heatmap of the predicted profiles for the best spacing and orientation,
	z-scored per sequence.
	"""

	# let's do the same for the predictions
	fig, ax = plt.subplots(1, 1, figsize = (8, 7))
	plt.set_loglevel('WARNING') 
	seaborn.heatmap(preds, cmap = "coolwarm", ax = ax, center = 0)

	# dashed lines to highlight edited regions
	ax.axvline(x = coords['start_A'] - SHIFT_REL_INPUT, color = "black", lw = 1, ls = "--")
	ax.axvline(x = coords['end_B'] - SHIFT_REL_INPUT, color = "black", lw = 1, ls = "--")

	plt.xlabel("Genomic position (1 kb)")
	plt.ylabel("Sequence")

	plt.xticks(
		ticks = [0, 100, 200, 300, 400, 500, 600, 700, 800, 900, 999],
		labels = [-500, -400, -300, -200, -100, 0, +100, +200, +300, +400, +500])

	# remove ticks
	plt.yticks([])
	plt.title(f"Predicted profile for {best_orientation} @ {best_spacing} bp (z-score per sequence)")
	plt.savefig(fig_path, dpi = 300);



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
	# compiled_modisco_h5_path = chrombpnet_dir + "02-compendium/modisco_compiled/modisco_compiled_filtered.h5"
	# modisco_compiled_unique = chrombpnet_dir + "03-syntax/01/motifs_compiled_unique.tsv"

	# output paths
	out = os.path.join(PROJ_OUT, "03-syntax/04b/in_silico_marginalization/")
	print(out)

	os.makedirs(out, exist_ok=True)


	# LOAD DATA --------------------------------------------------------------------
	logger.critical("Loading resources...")

	compo_to_test = pd.read_csv("04b-composites_to_test.tsv", sep="\t")

	# filter to where test_spacing == "Y"
	compo_to_test = compo_to_test[compo_to_test.test_spacing == "Y"]

	compo_to_test.head()

	# if debugging, use the provided motif
	if args.motif is not None:

		# check if the motif is in the list
		assert args.motif in compo_to_test.motif_name.values, f"Motif {args.motif} not found in the list of motifs to test."

		compo_to_test = compo_to_test[compo_to_test.motif_name == args.motif]
		logger.critical(f"DEBUG: Running with {compo_to_test}")


	# RUN COOPERATIVITY TESTS ------------------------------------------------------

	# for each composite motif, run in silico marginalziations at different spacings
	for row in compo_to_test.itertuples():

		datetimenow = datetime.datetime.now().strftime("%Y-%m-%d-%H:%M:%S")

		# setup --------------------------------------------------------------------
		# get info for the experiment
		comp_motif_name = row.motif_name
		comp_motif_name_safe = comp_motif_name.replace("/", ".")
		comp_type = row.category
		cluster = row.Cluster

		# make out directory and paths
		out_dir = os.path.join(out, comp_motif_name_safe)
		os.makedirs(out_dir, exist_ok=True)

		ch2 = logging.FileHandler(os.path.join(out_dir, comp_motif_name_safe + '-' + datetimenow + '.log'))
		ch2.setLevel(logging.WARNING)
		formatter = logging.Formatter('%(asctime)s %(levelname)-8s %(message)s')
		ch2.setFormatter(formatter)
		logger.addHandler(ch2)

		try:

			# check if done file exists
			if not args.overwrite and not args.downstream_only and not args.interpret and not args.interpret_only and os.path.exists(os.path.join(out_dir, ".done")):

				logger.critical(f"{comp_motif_name} ({comp_type}) already done, skipping")
				continue

			else:

				# write a message if done file exists
				if os.path.exists(os.path.join(out_dir, ".done")) and args.downstream_only:
					logger.critical(f"{comp_motif_name} is complete; overwriting downstream plots and summaries.")
				elif os.path.exists(os.path.join(out_dir, ".done")) and args.interpret_only:
					logger.critical(f"{comp_motif_name} is complete; running interpretation and downstream plots and summaries.")
				elif os.path.exists(os.path.join(out_dir, ".done")) and args.overwrite:
					logger.critical(f"{comp_motif_name} is complete; overwriting computations, plots, and summaries.")

				path_spacing_pkl = os.path.join(out_dir, "predictions.pkl")
				path_spacing_df = os.path.join(out_dir, "summarized_effects_per_spacing.tsv")
				path_contribs_all_pkl = os.path.join(out_dir, "contribs.pkl")
				path_contribs_mean_pkl = os.path.join(out_dir, "mean_contribs.pkl")
				path_solo_pkl = os.path.join(out_dir, "solo_predictions.pkl")
				path_solo_df = os.path.join(out_dir, "summarized_effects_solo.tsv")

				logger.critical(f"{comp_motif_name} ({comp_type})")

				# get sequences to test
				# pattern_name = get_pattern_name(comp_motif_name)
				# pattern_class = get_pattern_class(comp_motif_name)

				# trim CWM and get the high-affinity sequence for the composite
				# cwm_trimmed, _, _ = tanutils.utils.trim_cwm(cwm)
				# seqComp = tanutils.utils.get_high_affinity_seq(cwm_trimmed)

				seqA = row.test_seqA
				seqB = row.test_seqB

				if comp_type == "homocomposite":

					# test that seqA and seqB are the same otherwise raise error
					# assert seqA == seqB, "For homocomposite motifs, seqA and seqB must be the same."

					# A pair of the same motifs can have two different orientations.
					# - head-to-head
					# - tail-to-tail
					# - head to tail
					input_motifs = dict({
						"HT": [seqA, seqB],
						"HH": [seqA, reverse_complement(seqB)],
						"TT": [reverse_complement(seqA), reverse_complement(seqB)]
					})

				elif comp_type == "heterocomposite":

					# Two different motifs can have four orientations.
					# - A> B>: A, B, head-to-tail
					# - B> A>: B, A, head-to-tail
					# - A> B<: head-to-head
					# - A< B>: tail-to-tail
					input_motifs = dict({
						"A,B HT": [seqA, seqB],
						"B,A HT": [seqB, seqA],
						"HH": [seqA, reverse_complement(seqB)],
						"TT": [reverse_complement(seqA), seqB]
					})

				logger.critical(input_motifs)

				# load data and model

				# if the .done file exists, no more computation is needed;
				# don't load the model/backgrounds 
				if not args.downstream_only:

					logger.critical("\tLoading model...")
				
					# load models
					models = [
						ChromBPNetWrapper(BPNet.from_chrombpnet(filename=os.path.join(PROJ_IN, f"01-models/models/bias_Heart_c0_thresh0.4/{cluster}/fold_0/models/chrombpnet_nobias.h5"))),
						ChromBPNetWrapper(BPNet.from_chrombpnet(filename=os.path.join(PROJ_IN, f"01-models/models/bias_Heart_c0_thresh0.4/{cluster}/fold_1/models/chrombpnet_nobias.h5"))),
						ChromBPNetWrapper(BPNet.from_chrombpnet(filename=os.path.join(PROJ_IN, f"01-models/models/bias_Heart_c0_thresh0.4/{cluster}/fold_2/models/chrombpnet_nobias.h5"))),
						ChromBPNetWrapper(BPNet.from_chrombpnet(filename=os.path.join(PROJ_IN, f"01-models/models/bias_Heart_c0_thresh0.4/{cluster}/fold_3/models/chrombpnet_nobias.h5"))),
						ChromBPNetWrapper(BPNet.from_chrombpnet(filename=os.path.join(PROJ_IN, f"01-models/models/bias_Heart_c0_thresh0.4/{cluster}/fold_4/models/chrombpnet_nobias.h5")))
					]

					logger.critical("\tLoading backgrounds...")
					# load GC-matched background sequences
					nonpeaks_bed = os.path.join(PROJ_IN, f"00-inputs/gc_matched_negatives_subset/{cluster}__output_negatives_rand100.bed")

					X_nonpeaks = tangermeme.io.extract_loci(nonpeaks_bed, genome_fa, ignore = CHAR_IGNORE)
					logger.critical(X_nonpeaks.shape)
					logger.critical(X_nonpeaks[0])

					# ensure tensors are correctly typed and consistent before the function call
					X_nonpeaks = X_nonpeaks.float()
					logger.critical(X_nonpeaks.shape)

					# make a config dict
					# write a DF with config info
					config_df = pd.DataFrame({
						"motif_name": [comp_motif_name],
						"category": [comp_type],
						"cluster": [cluster],
						"n_spaces": [args.n_spaces],
						"seqA": [seqA],
						"seqB": [seqB]})
					
					# add the input motifs to the config df as strings
					for key, motif in input_motifs.items():
						config_df[key] = [f"{motif[0]}_{motif[1]}"]

					config_df.transpose()

				else:

					logger.critical("\tSkipping predictions.")
					config_df = pd.read_csv(os.path.join(out_dir, "config.tsv"), sep="\t")


				# ISM ----------------------------------------------------------------------

				# Spacing -------------------------------------------------------------------
				# check if the spacing predictions have already been computed
				if not args.overwrite and os.path.exists(path_spacing_pkl):

					logger.critical("\tLoading precomputed spacing ISM...")
					with open(path_spacing_pkl, "rb") as f:
						spacing_preds = pickle.load(f)
					
					spacing_dfs = pd.read_csv(path_spacing_df, sep="\t")

				else:

					logger.critical("\tRunning spacing ISM...")

					# store summary dataframes
					spacing_dfs = []

					# store predictions in a dict
					spacing_preds = {}

					for key, motifs in input_motifs.items():
						
						logger.critical(f"\t\t{key}")
						y_before, y_after, df = test_spacing(models, motifs, X_nonpeaks, key, n_spaces=args.n_spaces)

						spacing_preds[key] = {
							"y_before": y_before,
							"y_after": y_after
						}
						spacing_dfs.append(df)

					# reset index because when concatenating the dfs, they all
					# retain their original indices! and that will cause bugs
					# in getting the index of the max values when finding the 
					# best orientatino.
					spacing_dfs = pd.concat(spacing_dfs, ignore_index=True)

					# save the spacing predictions
					with open(path_spacing_pkl, "wb") as f:
						pickle.dump(spacing_preds, f)

					spacing_dfs.to_csv(path_spacing_df, sep="\t", index=False)

				# average the profiles
				# NOTE: this should be done the same way we consolidate profile values for generating the prediction tracks
				profile_preds_summarized = {}

				# for each motif in spacing preds dict, summarize the profile predictions
				for i, key in enumerate(input_motifs.keys()):

						y_before_profile_means, y_before_profile_sd = summarize_across_folds(mean_across_sequences(spacing_preds[key]["y_before"][0]))
						y_after_profile_means, y_after_profile_sd = summarize_across_folds(mean_across_sequences(spacing_preds[key]["y_after"][0]))

						profile_preds_summarized[key] = {
							"y_before_means": y_before_profile_means,
							"y_before_sd": y_before_profile_sd,
							"y_after_means": y_after_profile_means,
							"y_after_sd": y_after_profile_sd
						}
				
				# check the spacing df and find the spacing and orientation with the best effects fold change
				# which row maximizes the effect_mean?

				# DEBUG
				# import pdb
				# pdb.set_trace()

				max_idx = spacing_dfs.effect_mean.idxmax()
				max_row = spacing_dfs.iloc[max_idx]
				best_spacing = max_row.spacing
				best_orientation = max_row.key

				logger.critical(f"Best: {best_orientation} {best_spacing} bp, mean effect: {max_row['effect_mean']}")


				# solo effects -------------------------------------------------------------
				logger.critical("\tRunning solo motif ISM...")
				# Now, we do the ISM for composites and individual motifs alone
				# in silico marginalization for 3 sequences
				input_motifs_solo = {
					"A": seqA,
					"B": seqB
				}

				logger.critical(input_motifs_solo)

				# check if the spacing predictions have already been computed
				if not args.overwrite and os.path.exists(path_solo_pkl):

					logger.critical("\tLoading precomputed solo ISM...")
					with open(path_solo_pkl, "rb") as f:
						solo_preds = pickle.load(f)
					
					solo_dfs = pd.read_csv(path_solo_df, sep="\t")

				else:

					# store summary dataframes
					solo_dfs = []

					# store preds
					solo_preds = {}

					for key, motif in input_motifs_solo.items():

						logger.critical(f"\t\t{key}")
						y_before_list = []
						y_after_list = []

						for model in models:
							y_before, y_after = marginalize(model, X_nonpeaks, motif)
							y_before_list.append(y_before)
							y_after_list.append(y_after)

						# stack folds
						y_before = reshape_folds(y_before_list)
						y_after = reshape_folds(y_after_list)

						# calculate the mean counts across sequences, and then mean/sd across folds
						y_after_mean, y_after_sd = summarize_across_folds(mean_across_sequences(y_after[1]))

						# calculate the marginal effects in ln(counts) units
						y_effect = y_after[1] - y_before[1]

						# calculate mean across sequences, then the mean/sd across folds
						# thus the effect reported is a mean of differences
						y_effect_mean, y_effect_sd = summarize_across_folds(mean_across_sequences(y_effect))

						solo_preds[key] = {
							"y_before": y_before,
							"y_after": y_after
						}

						# create a dataframe
						df = pd.DataFrame({
							"key": [key],
							"counts_after_mean": y_after_mean.squeeze().tolist(),
							"counts_after_sd": y_after_sd.squeeze().tolist(),
							"effect_mean": y_effect_mean.squeeze().tolist(),
							"effect_sd": y_effect_sd.squeeze().tolist()
						})

						solo_dfs.append(df)

					# clean up
					del models
					torch.cuda.empty_cache()
					gc.collect()

					# combine the independent effects of A and B, handled differently for counts and profile
					# for profile, simply sum the effects
					y_before_profile_mult = solo_preds['A']['y_before'][0] 
					y_after_profile_mult = solo_preds['A']['y_after'][0] + solo_preds['B']['y_after'][0]

					# for counts, calculate the combined indepedent effects in the multiplicative model
					# counts log-additive effects
					y_before_counts_mult = solo_preds['A']['y_before'][1] + solo_preds['B']['y_before'][1]
					y_after_counts_mult = solo_preds['A']['y_after'][1] + solo_preds['B']['y_after'][1]
					
					# calculate the mean counts across sequences, and then mean/sd across folds
					y_before_counts_mean, y_before_counts_sd = summarize_across_folds(mean_across_sequences(y_before_counts_mult))
					y_after_counts_mean, y_after_counts_sd = summarize_across_folds(mean_across_sequences(y_after_counts_mult))

					# add to the dict
					solo_preds['A + B'] = {
						"y_before": [y_before_profile_mult, y_before_counts_mult],
						"y_after": [y_after_profile_mult, y_after_counts_mult]
					}

					# multiplicative models
					y_effect_mult = solo_preds['A']['y_after'][1] + solo_preds['B']['y_after'][1] - solo_preds['A']['y_before'][1] - solo_preds['B']['y_before'][1]

					# calculate mean across sequences, then the mean/sd across folds
					y_effect_mult_mean, y_effect_mult_sd = summarize_across_folds(mean_across_sequences(y_effect_mult))

					# sum A + B effects
					df_sum = pd.DataFrame({
						"key": ["A + B"],
						"counts_after_mean": y_after_counts_mean.squeeze().tolist(), 
						"counts_after_sd": y_after_counts_sd.squeeze().tolist(),
						"effect_mean": y_effect_mult_mean.squeeze().tolist(),
						"effect_sd": y_effect_mult_sd.squeeze().tolist()
					})

					solo_dfs.append(df_sum)
					solo_dfs = pd.concat(solo_dfs, ignore_index=True)
				
					# save the solo predictions
					with open(path_solo_pkl, "wb") as f:
						pickle.dump(solo_preds, f)

					solo_dfs.to_csv(path_solo_df, sep="\t", index=False)


				# average the profiles
				solo_profile_preds_summarized = {}

				# for each motif in spacing preds dict, summarize the profile predictions
				for i, key in enumerate(solo_preds.keys()):

						y_before_profile_means, y_before_profile_sd = summarize_across_folds(mean_across_sequences(solo_preds[key]["y_before"][0]))
						y_after_profile_means, y_after_profile_sd = summarize_across_folds(mean_across_sequences(solo_preds[key]["y_after"][0]))

						solo_profile_preds_summarized[key] = {
							"y_before_means": y_before_profile_means,
							"y_before_sd": y_before_profile_sd,
							"y_after_means": y_after_profile_means,
							"y_after_sd": y_after_profile_sd
						}

				logger.critical(solo_profile_preds_summarized['A']['y_before_means'].shape)

				# deep lift shap -----------------------------------------------------------
				# check if the deep lift shap values have already been computed
				if args.interpret or args.interpret_only:

					# don't overwrite, and the outputs exist
					if not args.overwrite_interpretation and os.path.exists(path_contribs_all_pkl) and os.path.exists(path_contribs_mean_pkl):

						logger.critical("\tLoading precomputed deep lift shap values...")
						with open(os.path.join(out_dir, path_contribs_all_pkl), "rb") as f:
							contributions = pickle.load(f)

						with open(os.path.join(out_dir, path_contribs_mean_pkl), "rb") as f:
							contributions_mean = pickle.load(f)

					else:

						logger.critical("Loading models with 5-fold wrapper for counts head...")

						# load model with five folds, using the CountsWrapper
						models_counts = [
							CountWrapper(BPNet.from_chrombpnet(filename=os.path.join(PROJ_IN, f"01-models/models/bias_Heart_c0_thresh0.4/{cluster}/fold_0/models/chrombpnet_nobias.h5"))),
							CountWrapper(BPNet.from_chrombpnet(filename=os.path.join(PROJ_IN, f"01-models/models/bias_Heart_c0_thresh0.4/{cluster}/fold_1/models/chrombpnet_nobias.h5"))),
							CountWrapper(BPNet.from_chrombpnet(filename=os.path.join(PROJ_IN, f"01-models/models/bias_Heart_c0_thresh0.4/{cluster}/fold_2/models/chrombpnet_nobias.h5"))),
							CountWrapper(BPNet.from_chrombpnet(filename=os.path.join(PROJ_IN, f"01-models/models/bias_Heart_c0_thresh0.4/{cluster}/fold_3/models/chrombpnet_nobias.h5"))),
							CountWrapper(BPNet.from_chrombpnet(filename=os.path.join(PROJ_IN, f"01-models/models/bias_Heart_c0_thresh0.4/{cluster}/fold_4/models/chrombpnet_nobias.h5")))
						]

						# make the edited sequences for the best orientation and spacing
						X_best = multisubstitute(X_nonpeaks, input_motifs[best_orientation], [best_spacing])
						X_best.shape

						# double check that the sequences are different
						logger.critical(characters(X_nonpeaks[0, :, 1040:1070]))
						logger.critical(characters(X_best[0, :, 1040:1070]))

						logger.critical("\tRunning deep lift shap...")
						contributions_list = []
						for model in models_counts:
							# run deep lift shap
							contributions_model = deep_lift_shap(model, X_best.float(), random_state=0, batch_size=16)
							contributions_list.append(contributions_model)
						contributions = torch.stack(contributions_list, dim=0)
						logger.critical(contributions.shape)

						# average across sequences and then folds
						contributions_mean, _ = summarize_across_folds(mean_across_sequences(contributions))
						contributions_mean.shape

						# save the spacing predictions
						with open(path_contribs_all_pkl, "wb") as f:
							pickle.dump(contributions, f)

						with open(path_contribs_mean_pkl, "wb") as f:
							pickle.dump(contributions_mean, f)

						# clean up
						del models_counts
						torch.cuda.empty_cache()
						gc.collect()


				# barplot for effect --------------------------------------------
				# create a figure with subplots in one row
				logger.critical("\tPlotting results...")
				logger.critical("\t\tplot_effect_spacing")
				plot_effect_spacing(spacing_dfs, input_motifs, fig_path = os.path.join(out_dir, f"effect_spacing.pdf"))
				plot_effect_spacing_lines(spacing_dfs, solo_dfs, comp_motif_name, fig_path = os.path.join(out_dir, f"effect_spacing_lines.pdf"))

				# solo effects -------------------------------------------------------------
				# make a small barplot of the solo effects
				logger.critical("\t\tplot_effect_solo")
				plot_effect_solo(solo_dfs, best_orientation, best_spacing, fig_path = os.path.join(out_dir, f"effect_solo.pdf"))


				# profile plots ------------------------------------------------------------
				logger.critical("\t\tplot_profiles_spacing")
				plot_profiles_spacing(profile_preds_summarized, input_motifs, fig_path = os.path.join(out_dir, f"profiles_per_spacing.pdf"))


				# solo profiles ------------------------------------------------------------
				logger.critical("\t\tplot_profiles_solo")
				plot_profiles_solo(solo_profile_preds_summarized, profile_preds_summarized, best_orientation, best_spacing, fig_path = os.path.join(out_dir, f"profiles_solo.pdf"))


				# calc effects across folds, per sequence -------------------------------------
				# effects of predicted counts across sequences for the best spacing/orientation, mean across folds
				
				# effects of predicted counts across sequences for the joint motifs at the joint spacing/orientation, mean across folds
				y_effect_joint = (
					spacing_preds[max_row.key]['y_after'][1][:, :, max_row.spacing, :] - spacing_preds[max_row.key]['y_before'][1][:, :, max_row.spacing, :]
					  ).mean(axis = 0)
				y_effect_joint.shape

				# effects of predicted counts across sequences for the summed effects, mean across folds
				y_effect_mult = (
					solo_preds['A']['y_after'][1] + solo_preds['B']['y_after'][1] - solo_preds['A']['y_before'][1] - solo_preds['B']['y_before'][1]
								).mean(axis = 0)

				y_effect_mult.shape

				# Z-scores & entropy -----------------------------------------------------------------
				# calculate the z-scores for the effects across spacings and
				# orientations. We want to test if effect at any spacing
				# significantly deviate from the mean
				effects = spacing_dfs.effect_mean.values

				# calculate Z scores
				Z = calc_effect_Z_score(effects)

				# calculate the entropy over the effects
				S = calc_effect_entropy(effects)

				# make a dataframe
				df_z = pd.DataFrame({
					"orientation": spacing_dfs.key,
					"spacing": spacing_dfs.spacing,
					"Z": Z})
				
				# save the Z-scores
				df_z.to_csv(os.path.join(out_dir, "Z_scores.tsv"), sep="\t", index=False)

				logger.critical("\t\tplot_effect_zscores")
				plot_effect_zscores(df_z, input_motifs, fig_path = os.path.join(out_dir, "Z_scores_dotplot.pdf"))

				# Wilcoxon signed-rank test -----------------------------------------------
				# H0: the effects are the same
				# H1: the difference between the joint and multiplicative effects is greater than 0
				results_dict = {
					"motif_name": comp_motif_name,
					# NOTE! seq1 and seq2 are the sequences as they're ordered in the best arrangement, not the input.
					# they may not be simply equal to seqA and seqB.
					"seq1": input_motifs[best_orientation][0],
					"seq2": input_motifs[best_orientation][1],
					"seqA_palindrome": is_palindrome(input_motifs[best_orientation][0]),
					"seqB_palindrome": is_palindrome(input_motifs[best_orientation][1]),
					"best_orientation": best_orientation,
					"best_spacing": best_spacing,
					"best_seq": f"{input_motifs[best_orientation][0]}{'N' * best_spacing}{input_motifs[best_orientation][1]}",
					"joint_effect": max_row.effect_mean,
					"joint_vs_sum_p_value": wilcoxon(y_effect_joint, y_effect_mult, alternative = "greater").pvalue,
					"max_Z": Z.max(),
					"entropy": S
				}

				# calculate center-to-center distance between the motifs
				dist_centers = len(input_motifs[best_orientation][0])/2 + best_spacing + len(input_motifs[best_orientation][1])/2
				results_dict["dist_centers"] = dist_centers

				logger.critical(results_dict)
				results = pd.DataFrame(results_dict)

				# plot attributions ---------------------------------------------------------

				# calculate coordinates for windows to plot
				coords = {
						# a window to zoom into (~ approx the center of the 2114 bp window)
						"zoom_start":  MIDPOINT - 50,
						"zoom_end":  MIDPOINT + 50
				}

				# windows for highlighting, relative to the midpoint of the 2114 bp window for which contributions are calculated
				if best_spacing % 2 == 0:
					coords["end_A"] = MIDPOINT - best_spacing/2 - 0.5
					coords["start_B"] = MIDPOINT + best_spacing/2 - 0.5
				elif best_spacing % 2 == 1:
					coords["end_A"] = MIDPOINT - best_spacing // 2 + 0.5
					coords["start_B"] = MIDPOINT + best_spacing // 2 + 1.5

				# windows for highlighting
				coords["start_A"] = coords["end_A"] - len(input_motifs[best_orientation][0])
				coords["end_B"] = coords["start_B"] + len(input_motifs[best_orientation][1])

				logger.critical(coords)

				# plots that require the contributions
				if args.interpret or args.interpret_only:

					# plot the mean attrib
					logger.critical("\t\tplot_mean_attributions")
					plot_mean_contributions(profile_preds_summarized, contributions_mean, coords, input_motifs, best_spacing, best_orientation,
							fig_path = os.path.join(out_dir, "contributions_mean.pdf"))

					# heatmap of attributions
					logger.critical("\t\tplot_contributions_heatmap")
					contribs_mean_across_folds, _ = summarize_across_folds(contributions)
					contribs_mean_across_folds.shape

					# sum across positions to get a single contribution score per position
					contribs_sum_mean_across_folds = contribs_mean_across_folds.sum(1)
					contribs_sum_mean_across_folds.shape

					plot_contributions_heatmap(contribs_sum_mean_across_folds, coords, best_orientation, best_spacing, fig_path = os.path.join(out_dir, "contributions_heatmap.pdf"))

				# predictions --------------------------------------------------------------
				spacing_preds_mean_across_folds, _ = summarize_across_folds(spacing_preds[best_orientation]['y_after'][0][:, :, best_spacing, :])
				spacing_preds_mean_across_folds = spacing_preds_mean_across_folds.squeeze(1)
				spacing_preds_mean_across_folds.shape

				# z-score the counts
				spacing_preds_mean_across_folds_z = []
				for i in range(100):
					spacing_preds_mean_across_folds_z.append((spacing_preds_mean_across_folds[i] - spacing_preds_mean_across_folds[i].mean(0)) / spacing_preds_mean_across_folds[i].std(0))

				# stack
				spacing_preds_mean_across_folds_z = torch.stack(spacing_preds_mean_across_folds_z)
				spacing_preds_mean_across_folds_z.shape

				plot_profile_heatmap_z(spacing_preds_mean_across_folds_z, coords, best_orientation, best_spacing, fig_path = os.path.join(out_dir, "profile_heatmap.pdf"))

				# plot the boxplots of the effects --------------------------------
				# make a dataframe with predicted log counts over the 100 examples
				df_pred_effect = pd.DataFrame({
					"example": np.arange(100),
					"Joint": y_effect_joint.squeeze(),
					"A + B": y_effect_mult.squeeze(),
					"A": (solo_preds['A']['y_after'][1] - solo_preds['A']['y_before'][1]).mean(axis = 0).squeeze(),
					"B": (solo_preds['B']['y_after'][1] - solo_preds['B']['y_before'][1]).mean(axis = 0).squeeze()
				})

				df_pred_effect.head()

				logger.critical("\t\tplot_effect_boxplot")
				plot_effect_boxplot(df_pred_effect, best_orientation, best_spacing, fig_path = os.path.join(out_dir, "effect_boxplot.pdf"))

				# close figures.
				plt.close("all")

				# save the outputs ---------------------------------------------------------

				# save the config if it doesn't exist
				if not os.path.exists(os.path.join(out_dir, "config.tsv")):
					config_df.to_csv(os.path.join(out_dir, "config.tsv"), sep="\t", index=False)
				
				# save the results
				results.to_csv(os.path.join(out_dir, "results.tsv"), sep="\t", index=False)

				# make an empty file ".done"
				open(os.path.join(out_dir, ".done"), 'a').close()

				logger.critical("\tDone with composite.")

		except Exception as e:
			logger.critical(f"Error processing {comp_motif_name}: {e}")
			logger.removeHandler(ch2)
			continue
			
		logger.removeHandler(ch2)

	logger.critical("Done.")



if __name__ == '__main__':
	args = parse_args()
	main(args)
