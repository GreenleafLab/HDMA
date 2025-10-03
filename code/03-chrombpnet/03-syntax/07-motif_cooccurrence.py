# Purpose: using ChromBPNet-derived motif hits, find motifs which are enriched
# for co-occurrence

import sys
import os
import argparse
import numpy as np
import pandas as pd

from collections import defaultdict
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests

# global variables ----------------------------------------------------

with open("../../ROOT_DIR.txt", 'r') as f:
	proj_dir = f.readline().strip()

with open("../../AK_PROJ_DIR.txt", 'r') as f:
	kundaje_dir = f.readline().strip()

# set up directories
chrombpnet_dir = os.path.join(proj_dir, "output/03-chrombpnet")
hits_dir = os.path.join(chrombpnet_dir, "02-compendium/hits_unified_motifs/reconciled_per_celltype_peaks")

# load annotated motifs
motif_annotation = pd.read_csv(os.path.join(chrombpnet_dir, "03-syntax/01/motifs_compiled_unique.tsv"), sep = "\t")
motif_annotation[["pattern_class", "idx_uniq", "motif_name", "pattern", "annotation", "category"]]

# output paths
out = os.path.join(proj_dir, "output/03-chrombpnet/03-syntax/07")
figout = os.path.join(proj_dir, "figures/03-chrombpnet/03-syntax/07")

NARROWPEAK_SCHEMA = ["chr", "start", "end", "1", "2", "3", "4", "5", "6", "summit"]



# arguments -----------------------------------------------------

def parse_args():
	
	parser = argparse.ArgumentParser()

	# space separated list of clusters
	parser.add_argument("--clusters", nargs='+', type=str, required=True,
		help="Space separated list of clusters.")
	args = parser.parse_args()
	print(args)

	return args



# helper functions -----------------------------------------------------

def load_data(cluster):
	"""
	Get peaks and hits for a given cluster
	"""

	print("@ loading peaks for cluster:", cluster)
	peaks_bed = os.path.join(chrombpnet_dir, f"00-inputs/chrombpnet_peaks/{cluster}__peaks_bpnet.narrowPeak.gz")

	# check that files exist
	assert os.path.exists(peaks_bed)

	X_peaks_coords = pd.read_csv(peaks_bed, sep='\t', header=None, names = NARROWPEAK_SCHEMA)

	X_peaks_coords.head()

	X_peaks_coords.shape
	max(X_peaks_coords.index)
	len(X_peaks_coords.index)

	print("@ loading hits for cluster:", cluster)

	# load hits
	hits = pd.read_csv(
		os.path.join(chrombpnet_dir,
				f"02-compendium/hits_unified_motifs/reconciled_per_celltype_peaks/{cluster}/counts_v0.23_a0.8_all/hits_unique.reconciled.annotated.tsv.gz"), sep="\t")

	# filter to only keep hits that are in the peaks
	hits = hits[hits['peak_id'].isin(X_peaks_coords.index)]
	max(hits['peak_id'])

	# add broad annotation
	hits = hits.merge(motif_annotation[['motif_name', 'category', 'annotation_broad']], on=['motif_name'], how='left')

	# filter to category == "base"
	hits = hits[hits['category'] == 'base']

	# print unique motifs
	print(np.unique(hits['motif_name']))

	return hits




def find_cooccurring_motifs(hits, motif_of_interest, correct = False):
	"""
	Find motifs that significantly co-occur with the motif of interest using Fisher's exact test.
	"""	

	# for each peak, get the set of unique motifs that occur in it
	peak_to_motifs = hits.groupby('peak_id')['annotation_broad'].unique().apply(list)

	# create a mapping from a motif to all the peaks it's in
	motif_to_peaks = defaultdict(set)
	for peak_id, motifs in peak_to_motifs.items():
		for motif in motifs:
			motif_to_peaks[motif].add(peak_id)

	all_motif_names = motif_to_peaks.keys()	
	all_peak_ids = hits['peak_id'].unique().tolist()
	
	results = []

	for partner_motif in all_motif_names:

		if partner_motif == motif_of_interest:
			continue

		# drawing peaks with instances of one motif, and estimating the probability that 
		# we see instances of another, modelled by the hypergeometric distribution:
		# N: total number of peaks w/ motif instances
		# n: number of peaks w/ motif of interest
		# K: total number of peaks w/ other motif
		# k: number of peaks with both motifs
		N = len(all_peak_ids)

		peaks_with_interest = motif_to_peaks.get(motif_of_interest, set())
		n = len(peaks_with_interest)

		peaks_with_partner = motif_to_peaks.get(partner_motif, set())
		K = len(peaks_with_partner)

		peaks_with_both = peaks_with_interest.intersection(peaks_with_partner)
		k = len(peaks_with_both)

		if k == 0 or n == 0 or K == 0:
			continue

		contingency_table = [[k, n - k], [K - k, N - n - K + k]]

		# 'greater' tests the hypothesis that the odds ratio is > 1 (enrichment)
		# the odds ratio is the ratio of:
		# - the odds of motif B occurring in peaks with motif A
		# - the odds of motif B occurring in peaks without motif A
		odds_ratio, p_value = fisher_exact(contingency_table, alternative='greater')

		results.append({
			'motif_A': motif_of_interest,
			'class_A': hits[hits['annotation_broad'] == motif_of_interest]['pattern_class'].iloc[0],
			'motif_B': partner_motif,
			'class_B': hits[hits['annotation_broad'] == partner_motif]['pattern_class'].iloc[0],
			'p_value': p_value,
			'odds_ratio': odds_ratio,
			'num_A': n,
			'num_B': K,
			'num_both': k
		})

	results_df = pd.DataFrame(results)

	if correct:
		p_values = results_df['p_value']
		reject, q_values, _, _ = multipletests(p_values, alpha=0.05, method='fdr_bh')
		results_df['q_value'] = q_values
		results_df['is_significant'] = reject
		co_motifs = results_df[results_df['is_significant']].sort_values(by='odds_ratio', ascending=False)
		print(f"Top co-occurring motifs for {motif_of_interest}:")
		return co_motifs
	else:
		return results_df


def compute_all_cooccurring_motifs(hits, cluster):
	"""
	Loop over all base motifs and run the co-occurrence enrichment analysis using Fisher's exact test.
	"""

	out_path = os.path.join(out, f"{cluster}__motif_cooccurrence_fisher_exact.tsv")


	
	all_motif_names = np.unique(hits['annotation_broad'])

	results_all = []

	# loop over each base motif (broad category)
	for motif in all_motif_names:
		print(motif)
		res = find_cooccurring_motifs(hits, motif, correct=False)
		results_all.append(res)

	results_all_df = pd.concat(results_all, ignore_index=True)

	# multiple testing correction
	p_values = results_all_df['p_value']
	reject, q_values, _, _ = multipletests(p_values, alpha=0.05, method='fdr_bh')
	results_all_df['q_value'] = q_values
	results_all_df['is_significant'] = reject

	# significant results
	cooc_motifs = results_all_df[results_all_df['is_significant']].sort_values(by='odds_ratio', ascending=False)

	cooc_motifs['log10q'] = -np.log10(cooc_motifs['q_value'] + 1e-300)  # add small value to avoid log(0)
	cooc_motifs['log2odds'] = np.log2(cooc_motifs['odds_ratio'])

	# assign cluster
	cooc_motifs['cluster'] = cluster

	# write to TSV
	print(f"@ writing results to {out_path}")
	cooc_motifs.to_csv(out_path, sep="\t", index=False)




def main(args):

	print(f"@ processing clusters: {args.clusters}")

	for cluster in args.clusters:
		
		print(f"@ cluster: {cluster}")

		out_path = os.path.join(out, f"{cluster}__motif_cooccurrence_fisher_exact.tsv")

		# check if the out file exists
		if os.path.exists(out_path):

			print(f"@ {out_path} exists, skipping...")

		else:

			hits_cluster = load_data(cluster)

			print("@ running co-occurrence")
			cooc_motifs = compute_all_cooccurring_motifs(hits_cluster, cluster)

			print("@ done.")



if __name__ == '__main__':
	args = parse_args()
	main(args)
