# Purpose: similar to `modisco report`, generate a report for the compiled/
# merged patterns with their CWM logos, and get matches to known motif
# database using tomtom.
#
# We use the functionality in the `report` module of tf-modisco-lite to 
# find the top matches to known motifs for each merged pattern. For the database of known 
# motifs, we use the set of motifs from Vierstra.
# The below section uses code from `modiscolite.report.report_motifs`, with a few customizations:
# - extract more info from tomtom report
# - add in additional info on the matched motifs using the Vierstra database/annotations

import sys
import h5py as h5
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os
import modiscolite.report
import types
from typing import List, Union
import tempfile
import shutil
import argparse



# parse args
def parse_args():
      parser = argparse.ArgumentParser(description='Get tomtom matches and make a custom report for a compiled modisco object.')
      parser.add_argument('--modisco-h5', type=str, help='Path to the compiled modisco h5 file.')
      parser.add_argument('--out-dir', type=str, help='Path to the output directory.')
      parser.add_argument('--meme-db', type=str, help='Path to the motifs database in meme format.')
      # parser.add_argument('--motif-anno', type=str, help='Path to the motif annotation file.')
      parser.add_argument('--trim-threshold', type=float, default=0.3, help='Probability threshold for trimming the PPM.')
      parser.add_argument('--trim-min-length', type=int, default=3, help='Minimum length of the trimmed PPM.')
      parser.add_argument('--n-matches', type=int, default=10, help='Number of top matches to fetch from tomtom.')
      parser.add_argument('--top-n-matches', type=int, default=3, help='Number of top matches whose logos to include in the report.')
      # parser.add_argument('--img-path-suffix', type=str, default='./', help='Path suffix for the image paths in the report.')
      # parser.add_argument('--is-writing-tomtom-matrix', type=bool, default=True, help='If True, write the tomtom matrix to a file.')
      parser.add_argument('--verbose', type=bool, default=False, help='If True, print more info.')

      args = parser.parse_args()
      print(args)

      # check if the paths are valid
      if not os.path.isfile(args.modisco_h5):
          raise ValueError(f'Invalid path to the modisco h5 file: {args.modisco_h5}')

      return args




def fetch_tomtom_matches(ppm, cwm, is_writing_tomtom_matrix, out_dir,
  	pattern_name, motifs_db, background=[0.25, 0.25, 0.25, 0.25],
  	tomtom_exec_path='tomtom', trim_threshold=0.3, trim_min_length=3):
    
    """Fetches top matches from a motifs database using TomTom.
    Args:
    	ppm: position probability matrix- numpy matrix of dimension (N,4)
    	cwm: contribution weight matrix- numpy matrix of dimension (N,4)
    	is_writing_tomtom_matrix: if True, write the tomtom matrix to a file
    	out_dir: directory for writing the TOMTOM file
    	pattern_name: the name of the pattern, to be used for writing to file
    	background: list with ACGT background probabilities
    	tomtom_exec_path: path to TomTom executable
    	motifs_db: path to motifs database in meme format
    	temp_dir: directory for storing temp files
    	trim_threshold: the ppm is trimmed from left till first position for which
    		probability for any base pair >= trim_threshold. Similarly from right.
    Returns:
    	list: a list of up to n results returned by tomtom, each entry is a
    		dictionary with keys 'Target ID', 'p-value', 'E-value', 'q-value'
    """
    
    _, fname = tempfile.mkstemp()
    _, tomtom_fname = tempfile.mkstemp()
    
    score = np.sum(np.abs(cwm), axis=1)
    trim_thresh = np.max(score) * trim_threshold  # Cut off anything less than 30% of max score
    pass_inds = np.where(score >= trim_thresh)[0]
    trimmed = ppm[np.min(pass_inds): np.max(pass_inds) + 1]
    
    # can be None of no base has prob>t
    if trimmed is None:
    	return []
    
    # trim and prepare meme file
    modiscolite.report.write_meme_file(trimmed, background, fname)
    
    if not shutil.which(tomtom_exec_path):
    	raise ValueError(f'`tomtom` executable could not be called globally or locally. Please install it and try again. You may install it using conda with `conda install -c bioconda meme`')
    
    # run tomtom
    cmd = '%s -no-ssc -oc . --verbosity 1 -text -min-overlap 5 -mi 1 -dist pearson -evalue -thresh 10.0 %s %s > %s' % (tomtom_exec_path, fname, motifs_db, tomtom_fname)
    os.system(cmd)
    
    # get e-value and query consensus as well as q-value
    print(tomtom_fname)
    tomtom_results = pd.read_csv(tomtom_fname, sep="\t", usecols=(1, 4, 5, 7))

    # TEMP
    # output_subdir = os.path.join(out_dir, "tomtom")
    # output_filepath = os.path.join(output_subdir, f"{pattern_name}.tomtom.tsv")
    # tomtom_results = pd.read_csv(output_filepath, sep="\t", usecols=(1, 4, 5, 7))

    os.system('rm ' + fname)
    if is_writing_tomtom_matrix:
    	output_subdir = os.path.join(out_dir, "tomtom")
    	os.makedirs(output_subdir, exist_ok=True)
    	output_filepath = os.path.join(output_subdir, f"{pattern_name}.tomtom.tsv")
    	os.system(f'mv {tomtom_fname} {output_filepath}')
    else:
    	os.system('rm ' + tomtom_fname)

    return tomtom_results



def main(args):

    # UNPACK ARGS -----------------------------------------------------
    modisco_h5 = args.modisco_h5
    out_dir = args.out_dir
    meme_db = args.meme_db
    trim_threshold = args.trim_threshold
    trim_min_length = args.trim_min_length
    n_matches = args.n_matches
    top_n_matches = args.top_n_matches
    verbose = args.verbose

    img_path_suffix = "./"
    is_writing_tomtom_matrix = True


    # LOAD MOTIFS --------------------------------------------------------
    motifs = modiscolite.report.read_meme(meme_db)
    len(motifs)

    # motifs_metadata = pd.read_csv(motif_anno, sep = "\t")
    # len(motifs_metadata)
    # motifs_metadata.head()


    # MAKE OUT DIRS ----------------------------------------------------
    modisco_logo_dir = os.path.join(out_dir, 'trimmed_logos')
    db_logo_dir = os.path.join(out_dir, 'db_logos')

    if not os.path.isdir(modisco_logo_dir):
        os.mkdir(modisco_logo_dir)

    if not os.path.isdir(db_logo_dir):
        os.mkdir(db_logo_dir)


    # MAKE LOGOS -----------------------------------------------------
    # load modisco obj and create fwd/rev trimmed logos per pattern
    pattern_groups = ['pos_patterns', 'neg_patterns']

    if verbose:
        print("@ making trimmed logos for modisco patterns...")

    modiscolite.report.create_modisco_logos(modisco_h5, modisco_logo_dir, trim_threshold, pattern_groups)

    # GET TOMTOM RESULTS ----------------------------------------------
    # initialize dataframe containing results
    results = {'pattern': [], 'modisco_cwm_fwd': [], 'modisco_cwm_rev': []}
    results

    with h5.File(modisco_h5, 'r') as modisco_results:
      for name in pattern_groups:
        if name not in modisco_results.keys():
          continue

        metacluster = modisco_results[name]
        for pattern_name, pattern in metacluster.items():
            
          results['pattern'].append(pattern_name)
          results['modisco_cwm_fwd'].append(os.path.join(img_path_suffix, 'trimmed_logos', f'{pattern_name}.cwm.fwd.png'))
          results['modisco_cwm_rev'].append(os.path.join(img_path_suffix, 'trimmed_logos', f'{pattern_name}.cwm.rev.png'))

    patterns_df = pd.DataFrame(results)
    reordered_columns = ['pattern', 'modisco_cwm_fwd', 'modisco_cwm_rev', 'query_consensus']
    patterns_df

    # get tomtom results
    tomtom_results = {}
    tomtom_results['query_consensus'] = []

    for i in range(n_matches):
        tomtom_results[f'match{i}'] = []
        # tomtom_results[f'TF{i}'] = []
        tomtom_results[f'qval{i}'] = []
        tomtom_results[f'e_val{i}'] = []
    
    if verbose:
        print("@ fetching tomtom matches...")

    with h5.File(modisco_h5, 'r') as modisco_results:
        for contribution_dir_name in pattern_groups:
            if contribution_dir_name not in modisco_results.keys():
                continue

            metacluster = modisco_results[contribution_dir_name]

            for pattern_name, pattern in metacluster.items():

                if verbose:
                    print("\t@ ", pattern_name)
                    
                ppm = np.array(pattern['sequence'][:])
                cwm = np.array(pattern['contrib_scores'][:]) 
            
                r = fetch_tomtom_matches(ppm, cwm, is_writing_tomtom_matrix, out_dir, pattern_name, meme_db)
                
                # get query target consensus sequence
                tomtom_results['query_consensus'].append(r.iloc[0]['Query_consensus'])

                # collect data up to n_matches rows -- each will be added as a new column in the eventual results df
                i = -1
                for i, (target, e_val, qval, _) in r.iloc[:n_matches].iterrows():
                    
                    tomtom_results[f'match{i}'].append(target)
                    tomtom_results[f'e_val{i}'].append(e_val)
                    tomtom_results[f'qval{i}'].append(qval)

                    # get the TF from the motif annotation
                    # tf = motifs_metadata[motifs_metadata['motif_id'] == target]['tf_name'].to_string(index = False).strip()
                    # tomtom_results[f'TF{i}'].append(tf)            

                # if there are fewer than n_matches rows of tomtom results, fill with None
                for j in range(i+1, n_matches):
                    tomtom_results[f'match{j}'].append(None)
                    # tomtom_results[f'TF{j}'].append(None)
                    tomtom_results[f'e_val{j}'].append(None)
                    tomtom_results[f'qval{j}'].append(None)		

    # convert to dataframe
    tomtom_results_df = pd.DataFrame(tomtom_results)
    patterns_df = pd.concat([patterns_df, tomtom_results_df], axis=1)

    if verbose:
        print("@ saving tomtom results...")    

    # MAKE KNOWN MOTIF LOGOS ------------------------------------------------------------------
    # make logos for motifs in database for the top N matches, and add paths to results data frame
    for i in range(top_n_matches):
        name = f'match{i}'
        
        logos = []

        for _, row in patterns_df.iterrows():

            if name in patterns_df.columns:
                if pd.isnull(row[name]):
                    logos.append("NA")
                else:
                    modiscolite.report.make_logo(row[name], db_logo_dir, motifs)
                    logos.append(f'./db_logos/{row[name]}.png')
            else:
                break

        patterns_df[f"{name}_logo"] = logos
        reordered_columns.extend([name, f'qval{i}', f'e_val{i}', f'{name}_logo'])

    # get column names for all other matches
    for j in range(i+1, n_matches):
        name = f'match{j}'
        reordered_columns.extend([name, f'qval{j}', f'e_val{j}'])

    print(reordered_columns)
    
    # reorder columns
    patterns_df = patterns_df[reordered_columns]

    # SAVE RESULTS ------------------------------------------------------
    patterns_df.to_csv(out_dir + "/modisco_compiled.tsv", sep="\t", index = False)

    print("@ done.")


if __name__ == '__main__':
    args = parse_args()
    main(args)
