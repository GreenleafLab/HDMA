# Purpose: split fragments (grouped per cell type) into pseudoreplicates for peak calling
# Adapted from Salil Deshpande. We assume that fragments are already split by cell type.

import os
import random
import sys

# specify allowed chromosomes
allowed_chros = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]

# parse the input file
filename = sys.argv[1]
cluster_sample = filename.split(".")[0] # remove extension
cluster = cluster_sample.split("__")[0] # get cluster
sample = cluster_sample.split("__")[-1] # get sample

print("--------------")
print(cluster_sample)
print(cluster)
print(sample)

# set up inputs and outputs
basedir = sys.argv[2]
frags_infile = f"{basedir}/fragments/{filename}"

out_pseudorep1_file = f"{basedir}/pseudorep1/{cluster}__{sample}.tsv"
out_pseudorep2_file = f"{basedir}/pseudorep2/{cluster}__{sample}.tsv"
out_pseudorepT_file = f"{basedir}/pseudorepT/{cluster}__{sample}.tsv"

# parse fragments files, splitting each one into pseudoreps
num_matches = 0
with open(frags_infile, 'r') as f_f_in, open(out_pseudorep1_file, 'w') as f_p1_out, open(out_pseudorep2_file, 'w') as f_p2_out, open(out_pseudorepT_file, 'w') as f_pT_out:
	for line in f_f_in:
		chro, start, end, barcode = line.strip().split("\t")
		start = int(start); end = int(end)
		if chro in allowed_chros:
			# Output pseudorepT
			f_pT_out.write(f"{chro}\t{start}\t{start+1}\t{barcode}\n")
			f_pT_out.write(f"{chro}\t{end-1}\t{end}\t{barcode}\n")
			# Output start to p1/p2
			if random.random() < 0.5:
				f_p1_out.write(f"{chro}\t{start}\t{start+1}\t{barcode}\n")
			else:
				f_p2_out.write(f"{chro}\t{start}\t{start+1}\t{barcode}\n")
			# Output end to p1/p2
			if random.random() < 0.5:
				f_p1_out.write(f"{chro}\t{end-1}\t{end}\t{barcode}\n")
			else:
				f_p2_out.write(f"{chro}\t{end-1}\t{end}\t{barcode}\n")

print(f"{sample} {cluster}")
