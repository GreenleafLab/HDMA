# Author: Selin Jessa
# 2024

import numpy as np
import matplotlib.pyplot as plt
import seaborn; seaborn.set_style('white')
import logomaker
import pandas as pd


INPUTLEN = 2114
OUTPUTLEN = 1000
SHIFT_REL_INPUT = np.int64((2114 - 1000) / 2)



def plot_motif(mat, type="CWM", title=None, ax=None, plot_IC=False, figsize=(6, 2), spines_visible = False):
	"""
	Plot a motif logo.

	Args:
		mat (np.array): The motif matrix, with shape (4, L) where L is the motif length.
	"""

	if ax is None:
		plt.figure(figsize=figsize)
		ax = plt.subplot(111)

	# do IC-scaling if type is PWM
	if type == "PPM" and plot_IC:
		mat_df = logomaker.transform_matrix(mat.T, from_type='probability', to_type='information')
		type = "PPM IC (bits)"
	else:
		mat_df = pd.DataFrame(mat.T, columns=["A", "C", "G", "T"])
	
	# tangermeme.plot.plot_logo(mat, ax=ax)
	logo = logomaker.Logo(mat_df, ax=ax)
	logo.style_spines(visible=spines_visible)
	
	ax.set_title(title)
	ax.set_ylabel(type)
	plt.tight_layout()
	
	if ax is None:
		plt.show()




def trim_cwm(cwm, trim_threshold=0.3, trim_min_length=3, flank=0):
    """
    Trim a CWM based on the contributions of the sequence to the pattern.
    This is adapted from Jacob Schreiber / tf-modiscolite
    (https://github.com/jmschrei/tfmodisco-lite/blob/3c6e38f3ad5df80c55bd4e8c7c2a531ee0a2b316/modiscolite/report.py#L721)
    and Ryan Zhao. The main difference between this and trim_ppm functions used is that here, we 
    trim based on the contributions directly, rather than on the probabilities. 

    Args:
        cwm (np.array): The CWM (contribution scores) of the sequence to the pattern.
        trim_threshold (float): The threshold for trimming; the minimum contribution score to keep a base.
        trim_min_length (int): The minimum length of the trimmed CWM.
    """
    score = np.sum(np.abs(cwm), axis=1)
    trim_thresh = np.max(score) * trim_threshold
    pass_inds = np.where(score >= trim_thresh)[0]

    # the start is the start of the passing idx - flank, or if the flank makes this extend into negative positions,
    # then the start is 0
    imp_start = np.max([np.min(pass_inds)-flank, 0])
    
    # the end is the end of the passing idx + flank, or if the flank makes this extend beyond the length of the seq,
    # then the end is the last position
    imp_end = np.min([np.max(pass_inds) + 1 + flank, cwm.shape[0]-1])
    trimmed = cwm[imp_start: imp_end]

    # can be None if no base has prob>t
    return trimmed, imp_start, imp_end



def get_high_affinity_seq(cwm):
	"""
	Function to get the highest affinity sequence from a CWM
	by taking the base with highest absolute contributions at each position.

	Args:
		cwm (np.array, shape L, 4): The CWM (contribution scores) of the sequence
	"""

	# get the highest affinity sequence
	high_affinity_seq = "".join([["A", "C", "G", "T"][np.argmax(np.abs(cwm[i]))] for i in range(cwm.shape[0])])

	return high_affinity_seq


HIT_DTYPE_DICT = {
	"seqnames": "str",
	"start": "int",
	"end": "int",
	"width": "int",
	"strand": "str",
	"start_untrimmed": "int",
	"end_untrimmed": "int",
	"motif_name": "str",
	"source": "str",
	"hit_coefficient": "float64",
	"hit_correlation": "float64",
	"hit_importance": "float64",
	"peak_name": "str",
	"peak_id": "int",
	"motif_name_unlabeled": "str",
	"pattern_class": "str",
	"distToGeneStart": "int",
	"nearestGene": "str",
	"peakType": "str",
	"distToTSS": "int",
	"nearestTSS": "str",
	"GC": "float64",
	"N": "int",
	"distToPeakSummit": "int"
}


def get_rel_hit_positions(hits, peaks, shift=SHIFT_REL_INPUT, flank=0):
	"""
	This function takes in a subset of a hits dataframe (which contain the absolute
	genomic coordinates of each hit), and the peak coordinatess ataframe, and gets the 
	positions of each hit relative to the input window (length 2114) and the output
	window (1000), to make it easy to do experiments on regions centered at hits.

	Returns
	-------
	np.ndarray
		The start positions of the hit relative to the input window (len 2114).
	np.ndarray
		The end positions of the hit relative to the input window (len 2114).
	np.ndarray
		The start positions of the hit relative to the output window (len 1000).
	np.ndarray
		The end positions of the hit relative to the output window (len 1000).

	"""

	print(hits.shape) 

	# get peak idx from hits
	peak_ids = hits['peak_id'].values
	peaks_subset = peaks.loc[peak_ids, :]

	# cast start, end columns in peaks df to int64
	peaks_subset['start'] = peaks_subset['start'].astype('int64')
	peaks_subset['end'] = peaks_subset['end'].astype('int64')

	# the relative start position of the motif within the peak
	hit_starts = hits['start'].values - peaks_subset['start'].values

	# the relative end position of the motif within the peak
	hit_ends = hits['end'].values - peaks_subset['start'].values

	# then get the hit start/end relative to the coordinates of the contribs, since they 
	# are length 2114 while the peaks and predictions are length 1000
	hit_starts2 = hit_starts + shift - flank
	hit_ends2 = hit_ends + shift + flank

	# make a dict
	hit_coords = {
		"start_input": hit_starts2,
		"end_input": hit_ends2,
		"start_output": hit_starts,
		"end_output": hit_ends
	}

	# return coords on the input len (2114), then coords on the output len (1000)
	return 	hit_coords