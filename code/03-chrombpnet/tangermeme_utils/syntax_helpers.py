
from .constants import *

import numpy as np
import pandas as pd

# plotting
import matplotlib.pyplot as plt
import seaborn; seaborn.set_style('white')
import matplotlib.ticker as ticker
import plotnine as pn
import logomaker
from tangermeme.plot import plot_logo

import bpnetlite
from bpnetlite import ChromBPNet
from bpnetlite import BPNet
from bpnetlite.bpnet import CountWrapper
from bpnetlite.attribute import deep_lift_shap

import torch

import tangermeme
from tangermeme.utils import random_one_hot 
from tangermeme.utils import one_hot_encode
from tangermeme.utils import characters
from tangermeme.utils import reverse_complement
from tangermeme.utils import _validate_input
from tangermeme.utils import _cast_as_tensor
from tangermeme.space import space
from tangermeme.ersatz import multisubstitute
from tangermeme.ersatz import substitute
from tangermeme.marginalize import marginalize

from scipy.stats import entropy

from plotnine import (
	ggplot, aes, geom_boxplot, geom_jitter, theme, xlab, ylab, ggtitle, element_blank, element_text, element_rect, scale_fill_manual, coord_cartesian, geom_hline
)

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

def plot_effect_spacing(spacing_dfs, input_motifs, fig_path, n_spaces):
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
		axes[i].set_xlim([0, n_spaces])
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

