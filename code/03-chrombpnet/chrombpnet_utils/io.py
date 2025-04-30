# Author: Selin Jessa
# 2024
# Purpose: Helper functions for loading model outputs, and outputs from downstream
# analysis (e.g. MoDISco, hit-calling)

import pandas as pd
import h5py
import numpy as np
from . import utils

NARROWPEAK_SCHEMA = ["chr", "start", "end", "1", "2", "3", "4", "5", "6", "summit"]

def load_peaks(peaks_path):

	peaks_df = pd.read_csv(peaks_path, header=None, sep='\t', names=NARROWPEAK_SCHEMA)
	return peaks_df



def load_hits(hits_path):

	dtype_dict = {'chr': 'str',
				'start': 'int',
				'end': 'int',
				'start_untrimmed': 'int',
				'end_untrimmed': 'int',
				'motif_name': 'str',
				'hit_coefficient': 'float64',
				'hit_correlation': 'float64',
				'strand': 'str',
				'peak_name': 'str',
				'peak_id': 'int'}

	hits_df = pd.read_csv(hits_path, sep='\t', dtype=dtype_dict)
	return hits_df




def load_reconciled_hits(hits_path):

	dtype_dict = {'seqnames': 'str',
				'start': 'int',
				'end': 'int',
				'start_untrimmed': 'int',
				'end_untrimmed': 'int',
				'motif_name': 'str',
				'source': 'str',
				'hit_coefficient': 'float64',
				'hit_correlation': 'float64',
				'strand': 'str',
				'peak_name': 'str',
				'peak_id': 'int',
				'motif_name_unlabeled' : 'str',
				'pattern_class' : 'str'}

	hits_df = pd.read_csv(hits_path, sep='\t', dtype=dtype_dict)
	return hits_df




def extract_from_bw(bw, chr, start, end):
	"""
	Extract data from a bigwig file in a given range.

	Args:
		bw: PyBigWig file handle
	"""
	
	scores = np.nan_to_num(bw.values(chr, start, end)).ravel()
	
	return scores




def reverse_complement(seq_ohe):
	"""
	Reverse complement a one-hot encoding of a DNA sequence, assuming the provided
	sequence is an array where the columns are A, C, G T
	"""
	
	# assuming ACGT,  the first reversal [::-1] reverses the sequence from 5' to 3',
	# and the second reverses elements in each sub array, so that
	# A <--> T since (1, 0, 0, 0) <--> (0, 0, 0, 1)
	# C <--> G since (0, 1, 0, 0) <--> (0, 0, 1, 0)
	return seq_ohe[::-1][:,::-1]




def get_hits_seqs(hits_df, genome, width=15, ohe=True, revcomp = True, revcomp_strand='-'):
	"""
	Fetches sequence from a given genome, given a df with chr, start, end, strand, strand.
	Specifically, we fetch the sequence centered at the midpoint between
	start, end, spanning a total length of 'width'. Optionally, reverse complement hits on
	one strand, and align seqlets so that the hit start is at the same position in both the
	forward and reverse-complemented sequences.

	Args:
		hits_df: pandas DataFrame with at least chr, start, end corresponding to motif hits
				 (could be from contribution scores/finemo, PWM-scanning, or other)
		genome: pyfaidx handle to reference genome Fasta file
	"""

	print("@ output width: " + str(width))
	
	vals = []

	for i, r in hits_df.iterrows():
		hit_width = r['end'] - r['start']

		# if reverse complementing, then if the hit_width is odd length, bump
		# center so that the nucleotides in the actual hit will be aligned
		# after the reverse complementing
		if hit_width % 2 == 1 and r['strand'] == revcomp_strand:
			center = 1 + r['start'] + hit_width//2
		else:
			center = r['start'] + hit_width//2
	   
		sequence = str(genome[r['chr']][(center - width//2):(center + width//2)])
		vals.append(sequence)

	if ohe:
		ohe = utils.dna_to_one_hot(vals)
		if revcomp:
			ohe = [reverse_complement(seq) if
				   hits_df['strand'].iloc[idx]==revcomp_strand else seq for idx, seq in enumerate(ohe)]
		return ohe
		
	else:
		return vals





def get_hits_conservation(hits_df, genome, conservation_bw, width=15, revcomp_strand = '-'):
	"""
	Fetches sequence from a given genome, given a df with chr, start, end, strand, strand.
	Specifically, we fetch the sequence centered at the midpoint between
	start, end, spanning a total length of 'width'. Optionally, reverse complement hits on
	one strand, and align seqlets so that the hit start is at the same position in both the
	forward and reverse-complemented sequences.

	Args:
		hits_df: pandas DataFrame with at least chr, start, end corresponding to motif hits
				 (could be from contribution scores/finemo, PWM-scanning, or other)
		genome: pyfaidx handle to reference genome Fasta file
		conservation_bw: str. Path to bigwig of conservation scores e.g. phyloP
	"""

	print("@ output width: " + str(width))

	import pyBigWig
	
	# open bigwig containing conservation
	cons_bigwig = pyBigWig.open(conservation_bw)
	
	conservation = []

	for i, r in hits_df.iterrows():
		
		hit_width = r['end'] - r['start']

		# if reverse complementing, then if the hit_width is odd length, bump
		# center so that the nucleotides in the actual hit will be aligned
		# after the reverse complementing
		if hit_width % 2 == 1 and r['strand'] == revcomp_strand:
			center = 1 + r['start'] + hit_width//2
		else:
			center = r['start'] + hit_width//2

		# coordinates of region to extract
		region_chr = r['chr']
		region_start = (center - width//2)
		region_end = (center + width//2)

		# get conservation
		cons_scores = extract_from_bw(cons_bigwig, region_chr, region_start, region_end)
		conservation.append(cons_scores)
		
	return conservation





def read_meme(filename, n_motifs=None):
	"""Read a MEME file and return a dictionary of PWMs.

	From tangermeme.io.read_meme: https://github.com/jmschrei/tangermeme/blob/main/tangermeme/io.py#L392
	This method takes in the filename of a MEME-formatted file to read in
	and returns a dictionary of the PWMs where the keys are the metadata
	line and the values are the PWMs.


	Parameters
	----------
	filename: str
		The filename of the MEME-formatted file to read in


	Returns
	-------
	motifs: dict
		A dictionary of the motifs in the MEME file.
	"""

	motifs = {}

	with open(filename, "r") as infile:
		motif, width, i = None, None, 0

		for line in infile:
			if motif is None:
				if line[:5] == 'MOTIF':
					motif = line.replace('MOTIF ', '').strip("\r\n")
				else:
					continue

			elif width is None:
				if line[:6] == 'letter':
					width = int(line.split()[5])
					pwm = np.zeros((width, 4))

			elif i < width:
				pwm[i] = list(map(float, line.strip("\r\n").split()))
				i += 1

			else:
				motifs[motif] = pwm.T
				motif, width, i = None, None, 0

				if n_motifs is not None and len(motifs) == n_motifs:
					break

	return motifs




# def load_preds(folds_paths):
	
#     y_pred_all = np.array([])
	
#     for fold in range(5):
	
#         print(f"@ fold {fold}")
	
#         # read in the predictions for each fold
#         preds_fold_path = folds_paths[fold]
#         preds_fold_h5 = h5py.File(preds_fold_path)
		
#         y_pred_fold = preds_fold_h5['predictions']['logcounts']
		
#         # concatenate the labels and predicted values for all folds
#         y_pred_all = np.concatenate([y_pred_all, y_pred_fold])
	
#     # check the shape of the concatenated arrays
#     print(y_pred_all.shape)
	
#     return y_pred_all


