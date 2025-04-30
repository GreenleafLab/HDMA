# Purpose: format peaks for chrombpnet. Recenter the peaks to the summits, and 
# adjust coordinates to 1,000 bp, maintaining a uniform region size for downstream analysis.
# Peaks will by 1,000 bp, centered on the peak summit.
# Adapted from Salil Deshpande.

import os
import sys

cluster = sys.argv[1]
peaks_dir = sys.argv[2]
out_dir = sys.argv[3]

print(f"----- {cluster} -----")
inpeak_file = f"{peaks_dir}/{cluster}__peaks_overlap_filtered.narrowPeak"
out_file = f"{out_dir}/{cluster}__peaks_bpnet.narrowPeak"
with open(inpeak_file, "r") as f_in, open(out_file, "w") as f_out:
	for line in f_in:
		line_split = line.strip().split("\t")
		chro, start, peak = line_split[0], line_split[1], line_split[9]
		midpoint = int(start) + int(peak)
		start_new = midpoint - 500; end_new = midpoint + 500
		assert(end_new - start_new == 1000)
		f_out.write(f"{chro}\t{start_new}\t{end_new}\t.\t.\t.\t.\t.\t.\t500\n")
