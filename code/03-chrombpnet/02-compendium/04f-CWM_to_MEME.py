# Purpose: convert CWMs to MEME format for easy loading in R, etc.
# Adapted from Jacob Schreiber, Ivy Evergreen in
# tf-modiscolite, https://github.com/jmschrei/tfmodisco-lite/blob/3c6e38f3ad5df80c55bd4e8c7c2a531ee0a2b316/modiscolite/io.py#L221C1-L273C29
# Major change is to write out positive as well as negative patterns.

import argparse
from collections import OrderedDict
import os
import textwrap

import h5py
import hdf5plugin

from typing import List, Literal, Union

import numpy as np
import scipy

from modiscolite import util
from modiscolite import meme_writer



if __name__=='__main__':
    
    parser = argparse.ArgumentParser(description='Convert CWMs to MEME format.')
    parser.add_argument("--f", "--filename", help='The name of the h5 file to read.', required=True, type=str)
    parser.add_argument("--d", "--datatype", help='The datatype to use for the MEME file.', required=True, type=str)
    parser.add_argument("--o", "--output", help='The name of the MEME file to write.', required=True, type=str)
    
    args = parser.parse_args()
    print(args)
    
    filename = args.f
    datatype = args.d
    output_filename = args.o

    alphabet = 'ACGT'
    writer = meme_writer.MEMEWriter(
        memesuite_version='5',
        alphabet=alphabet,
        background_frequencies='A 0.25 C 0.25 G 0.25 T 0.25'
    )

    with h5py.File(filename, 'r') as grp:

         for pattern_type in ['pos_patterns', 'neg_patterns']:
            if pattern_type in grp:
                for name, datasets in grp[pattern_type].items():
                    probability_matrix = datasets['contrib_scores'][:]
                    motif = meme_writer.MEMEWriterMotif(
                        name=name,
                        probability_matrix=probability_matrix,
                        source_sites=1,
                        alphabet=alphabet,
                        alphabet_length=4)

                    writer.add_motif(motif)
    
    writer.write(output_filename)

