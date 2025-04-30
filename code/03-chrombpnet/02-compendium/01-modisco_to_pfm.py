# Purpose: Given the PPMs computed by MoDISco (along with CWMs), we first trim the PPMs to regions
# of high probability and sufficient length. We then convert the PPMs to PFMs
# by multiplying the trimmed PPMs by the number of seqlets associated with the motif
# or pattern. The resulting PFMs are concatenated into one text file. If multiple
# modisco runs are provided, the motifs are concatenated.
# 
# Code adapted from:
# Surag Nair:
# https://github.com/kundajelab/surag-scripts/blob/b3babd9d2f2876d51fc86730edafdfac06778c17/modisco/modiscolite_to_pfm.py
# and Ryan Zhao.

import sys
import os
import h5py
import argparse
import numpy as np

# allow custom utils to be visible (https://stackoverflow.com/a/23891673)
sys.path.append(os.path.join(os.path.dirname(__file__), ".."))

# import custom helpers
from chrombpnet_utils.modisco import cwm_to_trimmed_ppm,ppm_to_gimme

if __name__=='__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-c", "--config", required=True, type=str, help="2-column TSV with the first column being the name of the model (e.g. cell type), and the second being the full path to the model's modisco h5 output file.")
    parser.add_argument("-o", "--output", required=True, type=str)
    parser.add_argument("-p", "--pattern_classes", default=['pos_patterns', 'neg_patterns'], type=str, help="Type of pattern classes to extract, either pos_patterns or neg_patterns or both.")
    parser.add_argument("-t", "--threshold", default=0.3, type=float)
    parser.add_argument("-ml", "--min_length", default=3, type=int, help="Min length of acceptable motif")
    args = parser.parse_args()
    
    # open the output file
    with open(args.output, 'w') as pfm_file:

        # for each line in the config file
        for line in open(args.config):
            
            # get tab-separated model name and modisco h5 file
            dataset, modisco_h5py = line.strip().split("\t")

            # open the modisco h5 file
            print(f"\t@ {dataset} loading {modisco_h5py}")
            f = h5py.File(modisco_h5py, 'r')
            
            # handle argument of length 1
            if isinstance(args.pattern_classes, str):
                args.pattern_classes = [args.pattern_classes]
                
            # print(args.pattern_classes)

            for pat_class in args.pattern_classes:

                if pat_class in f.keys():

                    for pattern in f[pat_class].keys():

                        # print(f"\t@ {pattern}")
                        
                        ppm = f[pat_class][pattern]['sequence']
                        cwm = f[pat_class][pattern]['contrib_scores']
                        trimmed_ppm = cwm_to_trimmed_ppm(ppm = ppm, cwm = cwm, trim_threshold=args.threshold, trim_min_length=args.min_length)
            
                        if trimmed_ppm is not None:
                            
                            # similar to the modisco reports, except prefixed by the name of the model / dataset
                            name = dataset + "__" + pat_class + "." + pattern
                            
                            # print(f"\t\t@ {name}")
                            pfm_file.write(ppm_to_gimme(name, trimmed_ppm))
                
            f.close()
            
