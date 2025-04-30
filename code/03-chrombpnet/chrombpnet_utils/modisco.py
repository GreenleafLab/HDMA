# Author: Selin Jessa, with some functions adapted from tf-modiscolite and Surag Nair
# (as indicated in the code.)
# 2024
# Purpose: helper functions for working with TF-MoDISco outputs


import numpy as np
import pandas as pd
import tqdm
import modiscolite
from . import io
from . import utils




def trim_ppm(ppm, t=0.45, min_length=3, flank=0):
    """
    Trim a PPM to include regions with high probability, which are also sufficiently long.
    Adapted from Surag Nair
    (https://github.com/kundajelab/surag-scripts/blob/b3babd9d2f2876d51fc86730edafdfac06778c17/modisco/modiscolite_to_pfm.py)

    Args:
        ppm (array): The input position-probability matrix
        t (float): The maximum probability threshold
        min_length (int): Minimum number of nucleotides required to keep a region
        flank (int): 

    """
    # trim matrix to first and last bp that have
    # p>=threshold 
    maxes = np.max(ppm,-1)
    # tuple of arrays representing row/col idxs with positions
    # above the threshold t
    maxes = np.where(maxes>=t) 

    # if no bases with prob>t or too small:
    if (len(maxes[0])==0) or (maxes[0][-1]+1-maxes[0][0]<min_length):
        return None

    # otherwise, returned the timmed PPM starting from the first base
    # that passes the threshold, optionally adding a certain amount
    # of flank on each side
    return ppm[max(maxes[0][0]-flank, 0):maxes[0][-1]+1+flank]



def trim_cwm(cwm, trim_threshold=0.3, trim_min_length=3, flank=0):
    """
    Trim a CWM based on the contributions of the sequence to the pattern.
    This is adapted from Jacob Schreiber / tf-modiscolite
    (https://github.com/jmschrei/tfmodisco-lite/blob/main/modiscolite/report.py#L72-L101)
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




def cwm_to_trimmed_ppm(ppm, cwm, trim_threshold=0.3, trim_min_length=3):
    """
    Trim a PPM based on the contributions of the sequence to the pattern.
    This is adapted from Jacob Schreiber / tf-modiscolite
    (https://github.com/jmschrei/tfmodisco-lite/blob/main/modiscolite/report.py#L72-L101)
    and Ryan Zhao. The main difference between this and trim_ppm functions used is that here, we 
    trim based on the contributions directly, rather than on the probabilities. 

    Args:
        ppm (np.array): The PPM (position probability matrix) to be trimmed.
        cwm (np.array): The CWM (contribution scores) of the sequence to the pattern.
        trim_threshold (float): The threshold for trimming; the minimum contribution score to keep a base.
        trim_min_length (int): The minimum length of the trimmed PPM.
    """
    score = np.sum(np.abs(cwm), axis=1)
    trim_thresh = np.max(score) * trim_threshold
    pass_inds = np.where(score >= trim_thresh)[0]
    trimmed = ppm[np.min(pass_inds): np.max(pass_inds) + 1]

    # can be None if no base has prob>t
    return trimmed








def ppm_to_gimme(name, ppm):
    """
    Adapted from Ryan Zhao
    Format a PPM to GIMME input format: https://gimmemotifs.readthedocs.io/en/master/reference.html#input-formats
    """
    result = f">{name}\n"
    for row in ppm:
        result += f"{row[0]:.4f}\t{row[1]:.4f}\t{row[2]:.4f}\t{row[3]:.4f}\n"
    return result






def get_seqlets_from_pattern(modisco_obj, pattern_name, peaks_path, pattern_type = "pos_patterns", window_size = 400):
    """
    For a given pattern in a given modisco run, extract the genomic coordinates of the
    seqlets that make up that pattern.
    Adapted from: tfmodiscolite.io.write_bed_from_h5 (https://github.com/jmschrei/tfmodisco-lite/blob/d98aeb17a9c79398ded3e21d0a30d004f39fdcd8/modiscolite/io.py#L276)

    Args:
        modisco_obj
        pattern
        pattern_type
        peaks
        window_size: int or None
    		The window size to use for the BED file. Default: 400, consistent with tf-modisolite.

    Returns:
        pandas DataFrame: data frame where each row is a seqlet, containing the columns:
            seqlet_chrom, seqlet_start, seqlet_end (genomic coordinates of the seqlets),
            strand, is_revcomp (indicating which strand the matching seqlet is found / if it was reverse-complemented for modisco),
            seqlet_id (indicating pattern type, pattern name, and the index of the seqlet),
            seqlet_coord (pretty version of coordinates, for pasting into genome browsers),
            peak_id (row index of the peak containing the seqlet, within the provided peaks file)

    """

    # extract the data associated with the pattern
    pattern_obj = modisco_obj[pattern_type][pattern_name]

    # load peaks
    with open(peaks_path, 'r') as peaks_file:
        peak_rows = peaks_file.read().splitlines()

    print(f"@ read in {str(len(peak_rows))} regions")
    
    seqlet_ids = []
    seqlet_chroms = []
    seqlet_genomic_starts = []
    seqlet_genomic_ends = []
    seqlet_strands = []
    is_revcomp = []
    peak_ids = []
    
    # process each seqlet within the pattern.
    for idx in tqdm.trange(pattern_obj['seqlets']['start'].shape[0]):
        
        seqlet_name = f'{pattern_name}.{idx}'
        # print(seqlet_name)
        
        row_num = pattern_obj['seqlets']['example_idx'][idx]
        peak_row = peak_rows[row_num].split('\t')
        chrom = peak_row[0]
        score = peak_row[4]
        
        # Seqlet starts and ends are offsets relative to the given
        # window, and the window's centers aligned with the peak's center.

        # Calculate the start and ends.
        absolute_peak_center = (int(peak_row[1]) + int(peak_row[2])) // 2

        window_center_offset = window_size // 2

        seqlet_start_offset = pattern_obj['seqlets']['start'][idx] + 1
        seqlet_end_offset = pattern_obj['seqlets']['end'][idx]
        absolute_seqlet_start = absolute_peak_center - window_center_offset + seqlet_start_offset
        absolute_seqlet_end = absolute_peak_center - window_center_offset + seqlet_end_offset

        strand_char = '-' if bool(pattern_obj['seqlets']['is_revcomp'][idx]) is True else '+'

        # append to our lists
        seqlet_ids.append(pattern_type + "." + seqlet_name)
        seqlet_chroms.append(chrom) 
        seqlet_genomic_starts.append(absolute_seqlet_start)
        seqlet_genomic_ends.append(absolute_seqlet_end)
        is_revcomp.append(bool(pattern_obj['seqlets']['is_revcomp'][idx]))
        seqlet_strands.append(strand_char)
        peak_ids.append(row_num)

    # make a nice dataframe
    seqlet_df = pd.DataFrame({
        "seqlet_chrom": seqlet_chroms,
        "seqlet_start": seqlet_genomic_starts,
        "seqlet_end": seqlet_genomic_ends,
        "strand": seqlet_strands,
        "is_revcomp": is_revcomp,
        "seqlet_id": seqlet_ids,
        "seqlet_coord": [chrom + ":" + str(start) + "-" + str(end) for chrom, start, end in zip(seqlet_chroms, seqlet_genomic_starts, seqlet_genomic_ends)],
        "peak_id": peak_ids
    })

    return seqlet_df






def extract_pattern(modisco_obj, pattern_name, peaks_path, contrib_scores, rev_comp = False, flank=0,
                    pattern_type = "pos_patterns", input_size = 2114, window_size = 400):
    """
    For a given pattern (_aka_ metacluster) in a given modisco run, extract all the relevant 
    characteristics and modelling outputs for that pattern and its seqlets.
    The function makes use of other helper functions that get each of those components separately.

    Args:
        modisco_obj
        pattern
        pattern_type
        peaks_path
        contrib_scores
        rev_comp: bool. Whether to reverse complement the CWM and all seqlet one-hot encodings
            and contribution scores (shaps), e.g. to match a known motifs.
        flank: int. How many bp to add on either side of the subset of the seqlets
            trimmed by importance. Passed to trim_cwm. Default: 0.
        input_size: int or None. The size (length in bp) of the input sequences.
            Default: 2114, consistent with chrombpnet.
        window_size: int or None. The size (length in bp) of the window surrounding the
            peak center that is considered for motif discovery.
            in modisco. Default: 400, consistent with tf-modisolite.
        
    """

    # extract the data associated with the pattern
    pattern_obj = modisco_obj[pattern_type][pattern_name]

    # get the contribution scores (CWM) for the pattern
    cwm = pattern_obj["contrib_scores"][()]

    # load peaks
    # TODO: this is done twice, once here and once within get_seqlets_from_pattern
    with open(peaks_path, 'r') as peaks_file:
        peaks = peaks_file.read().splitlines()
    
    # 1. trim the CWM by importance, to automatically get the start and end of the pattern
    trimmed_cwm, imp_start, imp_end = trim_cwm(cwm, flank=flank)

    # 2. get absolute coordinates of each seqlet  
    seqlet_df = get_seqlets_from_pattern(modisco_obj = modisco_obj,
                                         pattern_name = pattern_name,
                                         peaks_path = peaks_path,
                                         pattern_type = pattern_type,
                                         window_size = window_size)

    print(f"@ returning {seqlet_df.shape[0]} seqlets")
        
    # 3. get the one-hot encodings of each seqlet and the shap scores, trimmed to the pattern start/end,
    assert(len(peaks)==contrib_scores['shap']['seq'].shape[0])
    assert(input_size==contrib_scores['shap']['seq'].shape[-1])
    
    # region sliced out and used for modisco - generally this is 2114bp
    # i.e. start of the window is the midpoint of the sequence, minus half of the window size
    window_crop_start = input_size//2 - window_size//2
    window_crop_end = input_size//2 + window_size//2

    seqlet_one_hots = []
    seqlet_shaps = []
    
    # process each seqlet using the seqlet dataframe
    for index, row in seqlet_df.iterrows():

        # get the one-hot encoding, and then trim it using the trimmed pattern coordinates calculated above
        cur_ohe = pattern_obj['seqlets']['sequence'][index]# .transpose()
        cur_ohe_trimmed = cur_ohe[imp_start:imp_end, :]

        cur_shaps = pattern_obj['seqlets']['contrib_scores'][index]
        cur_shaps_trimmed = cur_shaps[imp_start:imp_end, :]

        # as per tfmodisco-lite, sequences/scores/etc have all been flipped in order to be in the
        # correct orientation, but sometimes we want to flip them to match a known motif
        if rev_comp:
            seqlet_one_hots.append(reverse_complement(cur_ohe_trimmed))
            seqlet_shaps.append(reverse_complement(cur_shaps_trimmed))
        else:
            seqlet_one_hots.append(cur_ohe_trimmed)
            seqlet_shaps.append(cur_shaps_trimmed)

    if rev_comp:
        cwm=reverse_complement(cwm)
        print("@ reverse complementing CWM, seqlet OHEs, and seqlet SHAPs")
    
    return pattern_obj, cwm, seqlet_df, imp_start, imp_end, seqlet_one_hots, seqlet_shaps




def extract_hit_data(modisco_obj, contribs_bw, conservation_bw, hits, genome, pattern_name, motif_name,
                     pattern_class = "pos_patterns", revcomp = True,
                     revcomp_strand = '-', trim = True, flank = 4, input_size = 2114, window_size = 400):
    """
    For a given motif in one cell type, extract all the hits, OHE sequences, and SHAP scores.

    Args:
        modisco_obj: h5py file handle for the Modisco h5 object
        contribs_bw: str. Path to contribution scores bigwig
        conservation_bw: str. Path to bigwig of conservation scores
        hits: str. Path to reconciled hits
        genome: pyfaidx Fasta file for the genome
        pattern_name: str. Name of pattern of interest, matching what's in the Modisco object.
        motif_name: str. Name of motif of interest, matching what's in the 'motif_name' col of the hits.
        pattern_class: str. one of "pos_patterns" or "neg_patterns"
        rev_comp: bool. Whether to reverse complement the CWM and all seqlet one-hot encodings
            and contribution scores (shaps), e.g. to match a known motifs.
        trim: bool. Whether to trim CWM.
        flank: int. How many bp to add on either side of the subset of the seqlets
            trimmed by importance. Passed to trim_cwm. Default: 4, matching the flank
            used to plot trimmed CWMs in Modisco:
            https://github.com/jmschrei/tfmodisco-lite/blob/d98aeb17a9c79398ded3e21d0a30d004f39fdcd8/modiscolite/report.py#L233
        input_size: int or None. The size (length in bp) of the input sequences.
            Default: 2114, consistent with chrombpnet.
        window_size: int or None. The size (length in bp) of the window surrounding the
            peak center that is considered for motif discovery.
            in modisco. Default: 400, consistent with tf-modisolite.
        
    """

    import pyBigWig
    
    # extract the data associated with the pattern
    pattern_obj = modisco_obj[pattern_class][pattern_name]

    # get the CWM for the pattern
    cwm = pattern_obj["contrib_scores"][()]

    # open bigwig containing contribution scores
    imp_bigwig = pyBigWig.open(contribs_bw)

    # open bigwig containing conservation
    # cons_bigwig = pyBigWig.open(conservation_bw)
    
    # 1. trim the CWM by importance, to automatically get the start and end of the pattern
    # potentially with some flanking nucleotides on either side
    if trim:
        trimmed_cwm, imp_start, imp_end = trim_cwm(cwm, flank=flank)
        width = imp_end - imp_start
        print(f"@ returning hits trimmed to width {width}")
    else:
        trimmed_cwm = cwm
        imp_start = 0
        imp_end = 29
        width = 30
        print(f"@ returning hits without trimming.")

    # 2. get hits and reset the row index
    hits_df = io.load_reconciled_hits(hits)
    hits_df_subset = hits_df.loc[hits_df["motif_name"] == motif_name].reset_index()

    print(f"@ returning {hits_df_subset.shape[0]} hits")
        
    # 3. in one pass through the hits, get (in a window around the hit)
    # - one-hot encodings of each hit
    # - SHAP scores
    # - conservation scores
    vals = []
    contribs = []
    # conservation = []
    
    for idx, row in hits_df_subset.iterrows():
        hit_width = row['end'] - row['start']

        # if reverse complementing, then if the hit_width is odd length, bump
        # center so that the nucleotides in the actual hit will be aligned
        # after the reverse complementing
        if hit_width % 2 == 1 and row['strand'] == revcomp_strand:
            center = 1 + row['start'] + hit_width//2
        else:
            center = row['start'] + hit_width//2

        # coordinates of region to extract
        region_chr = row['seqnames']
        region_start = (center - width//2)
        region_end = (center + width//2)

        # get sequences
        sequence = str(genome[region_chr][region_start:region_end])
        vals.append(sequence)

        # get SHAP scores
        scores = io.extract_from_bw(imp_bigwig, region_chr, region_start, region_end)
        contribs.append(scores)

        # get conservation
        # cons_scores = io.extract_from_bw(cons_bigwig, region_chr, region_start, region_end)
        # conservation.append(cons_scores)

    # close the bigwig file
    imp_bigwig.close()

    # 4. get OHE
    ohe = utils.dna_to_one_hot(vals)

    # reverse complement the sequence, and reverse the contributions, if a hit is on the
    # strand to rev-comp
    if revcomp:
        print("@ reverse complementing hit one-hot-encodings and SHAPs")
        
        ohe = [io.reverse_complement(seq) if
               hits_df_subset['strand'].iloc[idx]==revcomp_strand else seq for idx, seq in enumerate(ohe)]
        contribs = [scores[::-1] for scores in contribs]

    # so that the imp_start/imp_end for indexing the CWM match the length of 
    # the contents of ohe and contribs
    if width % 2 == 1:
        imp_end = imp_end - 1
    
    return cwm, imp_start, imp_end, hits_df_subset, ohe, contribs #, conservation







# get the relative coordinates to define each seqlet
def get_relative_coords(pattern_obj):
    """
    Args:
        patt_obj H5 group, corresponding to a pattern object or subpattern object (expected to contain "seqlets")
    """

    # for each seqlet in the provided pattern or subpattern, concatenate the example idx and rel start/end into a string.
    # This should be unique per seqlet.
    rel_coords = [str(idx) + ":" + str(start) + "-" + str(end) for idx,start,end in
     zip(pattern_obj['seqlets']['example_idx'][()],
         pattern_obj['seqlets']['start'][()],
         pattern_obj['seqlets']['end'][()])]

    return rel_coords






def get_subpattern_idx(pattern_obj):

    pattern_rel_coords = get_relative_coords(pattern_obj)
    
    # figure out which are the subpatterns in the object
    subpatterns = [key for idx,key in enumerate(pattern_obj.keys()) if key.startswith('subpattern')]

    # for each subpattern, get relative coords and find their indices
    all_idx_by_subpattern = []
    for sub_i in subpatterns:
        print(sub_i)
        sub_rel_coords = get_relative_coords(pattern_obj[sub_i])
        sub_idx = [i for i, x in enumerate(pattern_rel_coords) if x in sub_rel_coords]
        all_idx_by_subpattern = all_idx_by_subpattern + sub_idx
    
    print(len(all_idx_by_subpattern))
    return all_idx_by_subpattern
  
  
  


    
