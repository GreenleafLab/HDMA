# Purpose: generate a report after running finemo-gpu, with a unified motif set obtained
# from collapsing patterns from multiple modisco runs (each from different trained ChromBPNet models).

# Builds heavily on Austin Wang's finemo package, and the report generated there,
# - IO module: https://github.com/austintwang/finemo_gpu/blob/eb3ca5296fd61629a0b158225881623a9cbd92c0/src/finemo/data_io.py
# - evaluation module: https://github.com/austintwang/finemo_gpu/blob/eb3ca5296fd61629a0b158225881623a9cbd92c0/src/finemo/evaluation.py
# - reporting function: https://github.com/austintwang/finemo_gpu/blob/eb3ca5296fd61629a0b158225881623a9cbd92c0/src/finemo/main.py#L124
#
# Initially written for version 0.16 of the finemo package, updated for v0.23.

import sys
import os
import warnings
import argparse
import h5py
import polars as pl
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties
from tqdm import tqdm
import finemo.data_io as data_io
import finemo.evaluation as evaluation
from jinja2 import Template
import logging


def parse_args():
    parser = argparse.ArgumentParser(description="FiNeMo report")
    
    # general arguments
    parser.add_argument("--peaks", type=str, default=None, help="Path to the peak file (narrowPeak format).")
    parser.add_argument("--out-dir", type=str, default=None, help="Output directory for results.")
    parser.add_argument("--report", type=str, default=None, help="Path to HTML report template.")
    parser.add_argument("--modisco-region-width", type=int, default=400, help="The width of the region extracted around each peak summit.")
    parser.add_argument("--merged-motif-map", type=str, default=None,
                      help="Path to a tsv, expected to contain one row per component motif (i.e. motif discovered in each celltype). For each component motif, indicates the name of the merged pattern it belongs to.")
    parser.add_argument("--merged-motif-anno", type=str, default=None,
                      help="Path to the merged motif annotation file (tsv).")
    parser.add_argument("--alpha", type=float, default=0.08, help="alpha value for motif calling.")
    
    # cell-type specific arguments (hit calling)
    parser.add_argument("--celltype", type=str, default=None, help="Cell type in which hits are being called (corresponds to the name of the model, e.g. 'Eye_c11').")
    parser.add_argument("--celltype-modisco-h5", type=str, default=None,
                      help="Path to the cell-type specific modisco h5 file.")
    
    # unified hit calling arguments
    parser.add_argument("--unified-modisco-h5", type=str, default=None,
                      help="Path to the unified modisco h5 file, which was used for hitcalling.")
    parser.add_argument("--hits", type=str, default=None,
                      help="The `hits.tsv` output file from `finemo call-hits`.")
    parser.add_argument("--regions", type=str, default=None,
                      help="A .npz file of input sequences and contributions. Must be identical to the data used for hit calling and tfmodisco motif calling.") 
    parser.add_argument("--annotation-drop", type=str, default=None,
                      help="Report the number of hits to motifs with this annotation (in merged-motif-anno), e.g. 'low-quality'.") 
    args = parser.parse_args()
    print(args)

    # check that the files exit
    assert os.path.exists(args.merged_motif_map), f"File {args.merged_motif_map} not found."
    assert os.path.exists(args.merged_motif_anno), f"File {args.merged_motif_anno} not found."
    assert os.path.exists(args.celltype_modisco_h5), f"File {args.celltype_modisco_h5} not found."
    assert os.path.exists(args.unified_modisco_h5), f"File {args.unified_modisco_h5} not found."
    assert os.path.exists(args.hits), f"File {args.hits} not found."
    assert os.path.exists(args.regions), f"File {args.regions} not found."

    return args




def plot_hit_dist(occ_df, out_dir):
    """Plot the total hit distribution"""
    
    fig, ax = plt.subplots(figsize=(8, 4))
    
    unique, counts = np.unique(occ_df.get_column("total"), return_counts=True)
    freq = counts / counts.sum()
    num_bins = np.amax(unique) + 1
    x = np.arange(num_bins)
    y = np.zeros(num_bins)
    y[unique] = freq
    ax.bar(x, y)
    
    ax.set_xlabel("Motifs per peak")
    ax.set_ylabel("Frequency")
    
    output_path = os.path.join(out_dir, "total_hit_distribution.png")
    plt.savefig(output_path, dpi=150)
    plt.close(fig)



def get_unified_motif_map(merged_anno, merged_map, celltype):
    """Get dataframes that define the relationship between component patterns and merged/unified patterns."""

    map_unified = (
        merged_anno
        .join(
                merged_map, 
                on=["merged_pattern"],
                how="left"
            )
        # reformatting of column names for consistency
        .with_columns(
            pl.concat_str(
                [
                    pl.col("pattern_class"),
                    pl.col("merged_pattern")
                ],
                separator=".",
            ).alias("merged_pattern"))
        .with_columns(
            pl.concat_str(
                [
                    pl.col("pattern_class"),
                    pl.col("pattern")
                ],
                separator=".",
            ).alias("component_pattern_short"))
        .sort("n_seqlets", descending = True)
    )

    # filter to component patterns from this cell type
    map_per_celltype = map_unified.filter(pl.col("component_celltype").is_in([celltype]))
    
    # get only the anno per merged pattern
    map_unified_distinct = ( map_unified
                            .unique(subset=["merged_pattern", "annotation"], maintain_order=True)
                            .select(["merged_pattern", "annotation"]) )
    map_unified_distinct.shape

    return map_unified_distinct, map_per_celltype




def strip_class_from_motif_name(m):
    """Make a shortened version of the motif name for better aesthetics when plotting"""

    # split by period and keep elements from index 1 onwards
    parts = m.split(".")[1:]

    # rejoin the remaining parts with periods
    m_new = ".".join(parts)

    return m_new




def align_cwms(cwm1, cwm2):
    """
    Finds the alignment between two shifted or reverse-complemented PWMs using dot product.
    
    Args:
    cwm1: First CWM as a numpy array of size (4, 30).
    cwm2: Second CWM as a numpy array of size (4, 30).
    
    Returns:
    A tuple containing the shift amount and the maximum dot product.
    """
    
    max_dot = -float('inf')
    best_shift = 0
    
    # calculate dot product for all possible shifts
    shifts = range(-cwm2.shape[1], cwm2.shape[1])

    dot_prods = []
    for i, shift in enumerate(shifts):
        shifted_cwm2 = np.roll(cwm2, shift, axis=1)
        dot_prods.append(np.sum(cwm1 * shifted_cwm2))
        # print(f"shift: {shift}, dotproduct: {dot_prods[i]}")

    best_idx = np.flatnonzero(dot_prods == np.max(dot_prods))
    best_shifts = np.array(shifts)[best_idx]
    best_shift_idx = np.argmin(np.abs(best_shifts))
    
    return best_shifts[best_shift_idx], np.max(dot_prods)




def shift_hit_coords(seqlets_df, hits_filtered, map_per_celltype, motifs_discovered_from_celltype, modisco_obj, celltype_modisco_obj):
    """For each merged pattern which is also a component pattern, calculate its CWM alignment (shift and RC)
    to the corresponding component CWM, and then shift the coordinates accordingly so that hits coordinates
    (are more likely to) match seqlet coordinates.
    """

    # we only need to shift motifs when we have seqlets:
    
    # parition dataframes by motif
    hits_filtered_partitioned_by_motif_fromcelltype = (
        hits_filtered
        .filter(pl.col("motif_name").is_in(motifs_discovered_from_celltype))
        .collect()
        .partition_by("motif_name")
    )

    hits_filtered_partitioned_by_motif_new = (
        hits_filtered
        .filter(pl.col("motif_name").is_in(motifs_discovered_from_celltype).not_())
        .collect()
        .partition_by("motif_name")
    )

    len(hits_filtered_partitioned_by_motif_fromcelltype)
    len(hits_filtered_partitioned_by_motif_new)

    # initialize the new, filtered motifs with the ones we don't need to shift.
    hits_filtered_partitioned_by_motif_shifted = hits_filtered_partitioned_by_motif_new

    for df in hits_filtered_partitioned_by_motif_fromcelltype:

        # get pattern info
        df_m = df
        m = df_m.select(pl.col("motif_name")).to_series().to_list().pop(0)
        m_matched = map_per_celltype.filter(pl.col("merged_pattern") == m)
        pattern_class = m_matched.select(pl.col("pattern_class")).to_series().to_list().pop(0)
        pattern = m_matched.select(pl.col("pattern")).to_series().to_list().pop(0)

        print(f"\tmotif: {m}, class: {pattern_class}, pattern: {pattern}")

        # get cwms
        cwm1 = modisco_obj[pattern_class][strip_class_from_motif_name(m)]['contrib_scores'][()].transpose()
        cwm2 = celltype_modisco_obj[pattern_class][pattern]['contrib_scores'][()].transpose()

        # check if revcomp, find shift
        shift_fc, best_dot_fc = align_cwms(cwm1, cwm2)
        shift_rc, best_dot_rc = align_cwms(cwm1, cwm2[::-1,::-1])

        if best_dot_rc > best_dot_fc:
            is_revcomp = True
            shift = shift_rc
        else:
            is_revcomp = False
            shift = shift_fc

        print(f"\tis_revcomp: {is_revcomp}, shift: {shift}")

        # now we have to do some more careful shifting of coordinates, and swapping the is_revcomp boolean:
        if shift == 0:
            if is_revcomp:
                df_m_shifted = df_m.with_columns(is_revcomp = ~pl.col("is_revcomp"))
            else:
                df_m_shifted = df_m
        elif shift > 0 and is_revcomp:
            df_m_shifted = (
                df_m.with_columns(
                    is_revcomp = ~pl.col("is_revcomp"),
                    start_untrimmed = (pl.when(pl.col("is_revcomp") == True)
                                       .then(pl.col("start_untrimmed") + (-shift))
                                       .otherwise(pl.col("start_untrimmed") + shift)
                ),
                    end_untrimmed = (pl.when(pl.col("is_revcomp") == True)
                                       .then(pl.col("end_untrimmed") + (-shift))
                                       .otherwise(pl.col("end_untrimmed") + shift)
                )
                )
            )
        else:
            df_m_shifted = df_m.with_columns(
                (pl.col("start_untrimmed") + shift),
                (pl.col("end_untrimmed") + shift)
                )

        hits_filtered_partitioned_by_motif_shifted.append(df_m_shifted)

    print(len(hits_filtered_partitioned_by_motif_shifted))
    hits_filtered_shifted = pl.concat(hits_filtered_partitioned_by_motif_shifted, how = "vertical_relaxed")
    print(hits_filtered_shifted.shape)

    return hits_filtered_shifted



def get_recall(motif_names_called, seqlets_df, hits_unique, hits_filtered_shifted, modisco_obj, motif_width, map_unified_distinct, regions):
    """Calculate statistics to evaluate hitcalling, including recall with respect to seqlets, and get CWMs."""

    seqlets_df_collected = seqlets_df.collect()
    
    print("@ seqlets df:")
    print(seqlets_df_collected.shape)
    print(seqlets_df_collected.head())

    print("@ hits df:")
    print(hits_filtered_shifted.shape)
    print(hits_filtered_shifted.head())

    # convert hits_filtered_shifted start_untrimmed, end_untrimmed to u32
    hits_filtered_shifted = hits_filtered_shifted.with_columns(
        pl.col("start_untrimmed").cast(pl.UInt32),
        pl.col("end_untrimmed").cast(pl.UInt32)
    )

    # get the hit/seqlet overlap, now using the shifted coordinates:
    overlaps_df_shifted = (
            # TODO: bug here?
            # polars.exceptions.ComputeError: datatypes of join keys don't match -
            # `start_untrimmed`: i64 on left does not match `start_untrimmed`: u32 on right
            hits_filtered_shifted.join(
                seqlets_df_collected, 
                on=["chr", "start_untrimmed", "is_revcomp", "motif_name"],
                how="inner",
                # validate='1:1'
            )
        )
    print(overlaps_df_shifted.shape)
    print(len(set(overlaps_df_shifted.select(pl.col("motif_name")).to_series().to_list())))

    overlaps_by_motif = overlaps_df_shifted.partition_by("motif_name", as_dict = True)

    # regions that were seqlets, but not hits
    seqlets_only_df = (
        seqlets_df_collected.join(
            hits_filtered_shifted, 
            on=["chr", "start_untrimmed", "is_revcomp", "motif_name"],
            how="anti",
            # validate='1:1'
        )
    )

    # regions that were hits but not seqlets
    hits_only_filtered_df = (
        hits_filtered_shifted.join(
            seqlets_df_collected, 
            on=["chr", "start_untrimmed", "is_revcomp", "motif_name"],
            how="anti",
            # validate='1:1'
        )
    )

    # split dataframes by motif names
    hits_by_motif = hits_unique.collect().partition_by("motif_name", as_dict=True)
    hits_fitered_by_motif = hits_filtered_shifted.partition_by("motif_name", as_dict=True)
    seqlets_by_motif = seqlets_df.collect().partition_by("motif_name", as_dict=True)
    overlaps_by_motif = overlaps_df_shifted.partition_by("motif_name", as_dict = True)
    seqlets_only_by_motif = seqlets_only_df.partition_by("motif_name", as_dict=True)
    hits_only_filtered_by_motif = hits_only_filtered_df.partition_by("motif_name", as_dict=True)

    print(len(motif_names_called))
    print(len(hits_by_motif))
    print(len(hits_fitered_by_motif))
    print(len(overlaps_by_motif))

    recall_data = {}
    cwms = {}
    dummy_df = overlaps_df_shifted.clear()

    # get the stats for each motif
    for m in tqdm(motif_names_called):
    
        # print(m)
    
        if m.startswith("pos"):
            pattern_class = "pos_patterns"
        elif m.startswith("neg"):
            pattern_class = "neg_patterns"
    
        # make a new data frame using the cols in dummy_df, for the calls associated with motif m
        hits = hits_by_motif.get(m, dummy_df)
        hits_filtered = hits_fitered_by_motif.get(m, dummy_df)
        seqlets = seqlets_by_motif.get(m, dummy_df)
        overlaps = overlaps_by_motif.get(m, dummy_df)
        seqlets_only = seqlets_only_by_motif.get(m, dummy_df)
        hits_only_filtered = hits_only_filtered_by_motif.get(m, dummy_df)
    
        recall_data[m] = {
            "seqlet_recall": 0 if seqlets.is_empty() else np.float64(overlaps.height) / seqlets.height,
            "num_hits_total": hits.height,
            "num_hits_restricted": hits_filtered.height,
            "num_seqlets": seqlets.height,
            "num_overlaps": overlaps.height,
            "num_seqlets_only": seqlets_only.height,
            "num_hits_restricted_only": None if seqlets.is_empty() else hits_only_filtered.height
        }

        # get CWMs
        cwms[m] = {
                # this is the compiled modisco obj
                "modisco_fc": modisco_obj[pattern_class][strip_class_from_motif_name(m)]['contrib_scores'][()].transpose(),
                "hits_fc": evaluation.get_cwms(regions, hits, motif_width)
            }
        
        cwms[m]["hits_rc"] = cwms[m]["hits_fc"][::-1,::-1]
        cwms[m]["modisco_rc"] = cwms[m]["modisco_fc"][::-1,::-1]
        
        # deviating from finemo, we want to here calculate the correlation
        # between the hits-derived CWM and the modisco CWM. They should be close.
        hits_cwm = cwms[m]["hits_fc"]
        modisco_cwm = cwms[m]["modisco_fc"]
        hnorm = np.sqrt((hits_cwm**2).sum())
        mnorm = np.sqrt((modisco_cwm**2).sum())
        cwm_cor = (hits_cwm * modisco_cwm).sum() / (hnorm * mnorm)
    
        recall_data[m]["hits_modisco_cwm_correlation"] = cwm_cor
    
        if not seqlets.is_empty():
    
            # add additional CWMS
            cwms[m]["seqlets_fc"] = evaluation.get_cwms(regions, seqlets, motif_width)
            cwms[m]["seqlets_only"] = evaluation.get_cwms(regions, seqlets_only, motif_width)
            cwms[m]["hits_restricted_only"] = evaluation.get_cwms(regions, hits_only_filtered, motif_width)
    
            # now we calculate correlated btwn restricted hits and seqlets
            hits_only_cwm = cwms[m]["hits_restricted_only"]
            seqlets_cwm = cwms[m]["seqlets_fc"]
            hnorm = np.sqrt((hits_only_cwm**2).sum())
            snorm = np.sqrt((seqlets_cwm**2).sum())
            cwm_cor = (hits_only_cwm * seqlets_cwm).sum() / (hnorm * snorm)
        
            recall_data[m]["hits_restricted_seqlets_cwm_correlation"] = cwm_cor

    print(len(recall_data))
    print(len(cwms))

    # convert to dataframe
    records = [{"motif_name": k} | v for k, v in recall_data.items()]
    print(map_unified_distinct.shape)
    recall_df = (
        pl.from_dicts(records).sort(pl.col("num_hits_total"), descending = True)
        .join(
                map_unified_distinct.rename({"merged_pattern": "motif_name"}), # rename the column to join on 
                on=["motif_name"],
                how="left")
    )
    
    print(recall_df.shape)
    print(recall_df.head())

    return recall_data, recall_df, cwms



def plot_hits_vs_restricted(recall_df, motifs_discovered_from_celltype, out_dir):
    """Plot a scatterplot of the number of hits vs restricted hits"""

    # plot the total number of *restricted* hits per motif
    fig, ax = plt.subplots(figsize=(6, 6))
    
    x = recall_df.select(pl.col("num_hits_total")).to_series().to_list()
    y = recall_df.select(pl.col("num_hits_restricted")).to_series().to_list()
    
    plt.scatter(x, y)
    ax.set_xlabel("Total hits")
    ax.set_ylabel("Restricted hits")
    
    lim = max(np.amax(x), np.amax(y))
    ax.axline((0, 0), (lim, lim), color="0.3", linewidth=0.7, linestyle=(0, (5, 5)))
    
    output_path = os.path.join(out_dir, "hits_vs_restricted_hits.png")
    plt.savefig(output_path, dpi=150)
    plt.close(fig)




def plot_hits_per_motif(recall_df, motifs_discovered_from_celltype, out_dir):
    """Plot the number of hits for each motif, on a log10 scale"""

    # plot the total number of hits per motif
    fig, ax = plt.subplots(figsize=(25, 7))
    
    motif_labels = recall_df.select(pl.col("motif_name")).to_series().to_list()
    colors = [("green" if motif in motifs_discovered_from_celltype else "gray") for motif in motif_labels]

    x = (recall_df
         .with_row_count()
         .with_columns(
            pl.concat_str(
                [
                    # get row index so that we can get unique names
                    pl.col("row_nr"),
                    pl.col("annotation")
                ],
                separator=": ",
            ).alias("new_annotation"))
         .select(pl.col("new_annotation"))
         .to_series().to_list()
        )

    # DEBUG:
    # df_missing = ( recall_df
    #               .filter(
    #                   pl.any_horizontal(pl.col("annotation").is_null())
    #               )
    #              )
    # print(df_missing)

    print(x)
    y=recall_df.select(pl.col("num_hits_total")).to_series().to_list()
    
    plt.bar(x = x,
            height = y,
            color = colors)
    ax.set_xlabel("Motif")
    ax.set_ylabel("# hits")
    ax.set_yscale('log')
    ax.tick_params(axis='x', labelrotation=90, labelsize=5)
    plt.tight_layout() # https://matplotlib.org/stable/users/explain/axes/tight_layout_guide.html
    
    output_path = os.path.join(out_dir, "hits_per_motif.png")
    plt.savefig(output_path, dpi=150)
    plt.close(fig)



def plot_hits_vs_seqlets(recall_data, motifs_discovered_from_celltype, out_dir):
    """Plot a scatterplot of hit counts vs seqlet counts"""

    x = []
    y = []
    m = []
    
    for k, v in recall_data.items():
    
        if k in motifs_discovered_from_celltype:    
            x.append(v["num_hits_total"])
            y.append(v["num_seqlets"])
            m.append(k)
    
    lim = max(np.amax(x), np.amax(y))
    
    fig, ax = plt.subplots(figsize=(8,8))
    ax.axline((0, 0), (lim, lim), color="0.3", linewidth=0.7, linestyle=(0, (5, 5)))
    ax.scatter(x, y, s=7)
    
    for i, txt in enumerate(m):
        ax.annotate(i, (x[i], y[i]), fontsize=12, weight="bold")
    
    ax.set_yscale('log')
    ax.set_xscale('log')
    
    ax.set_xlabel("Hits per motif")
    ax.set_ylabel("Seqlets per motif")
    
    output_path = os.path.join(out_dir, "hits_vs_seqlet_counts.png")
    plt.savefig(output_path, dpi=300)
    plt.close(fig)



def plot_hit_distributions(occ_df, motif_names, out_dir):
    """Plot distribution of hits per peak, for each of the motifs."""
    
    motifs_dir = os.path.join(out_dir, "motif_hit_distributions")
    os.makedirs(motifs_dir, exist_ok=True)

    for m in motif_names:
        fig, ax = plt.subplots(figsize=(6, 2))

        unique, counts = np.unique(occ_df.get_column(m), return_counts=True)
        freq = counts / counts.sum()
        num_bins = np.amax(unique) + 1
        x = np.arange(num_bins)
        y = np.zeros(num_bins)
        y[unique] = freq
        ax.bar(x, y)

        output_path = os.path.join(motifs_dir, f"{m}.png")
        plt.savefig(output_path, dpi=300)

        plt.close(fig)




def calc_cooccurrence(num_peaks, num_motifs, occ_df, motif_names_called, top_k_motif_names, k, out_dir):
    """Calculate motif co-occurrence, and then plot a matrix to visualize it for the top k motifs."""

    # # create empty occurrence matrix, then populate
    # occ_mat = np.zeros((num_peaks, num_motifs), dtype=np.int16)
    # for i, m in enumerate(motif_names_called):
    #     occ_mat[:,i] = occ_df.get_column(m).to_numpy()
        
    # # matrix multiplication to get cooccurrence [slow step]
    # occ_bin = (occ_mat > 0).astype(np.int32)
    # coocc = occ_bin.T @ occ_bin
    # occ_mat.shape
    # coocc.shape

    occ_path = os.path.join(out_dir, "motif_occurrences.tsv")
    data_io.write_occ_df(occ_df, occ_path)

    # repeat for top k motifs only.
    occ_mat_k = np.zeros((num_peaks, k), dtype=np.int16)
    for i, m in enumerate(top_k_motif_names):
        occ_mat_k[:,i] = occ_df.get_column(m).to_numpy()
        
    # matrix multiplication to get cooccurrence
    occ_bin_k = (occ_mat_k > 0).astype(np.int32)
    coocc_k = occ_bin_k.T @ occ_bin_k

    # plot the co-occurrence matrix
    cov_norm = 1 / np.sqrt(np.diag(coocc_k))
    matrix = coocc_k * cov_norm[:,None] * cov_norm[None,:]
    
    fig, ax = plt.subplots(figsize=(15, 15))
    
    # Plot the heatmap
    ax.imshow(matrix, interpolation="nearest", aspect="auto", cmap="Greens")
    
    # Set axes on heatmap
    ax.set_yticks(np.arange(len(top_k_motif_names)))
    ax.set_yticklabels(top_k_motif_names, size=6)
    ax.set_xticks(np.arange(len(top_k_motif_names)))
    ax.set_xticklabels(top_k_motif_names, rotation=90, size=6)
    ax.set_xlabel("Motif i")
    ax.set_ylabel("Motif j")
    
    fig.tight_layout()
    
    output_path = os.path.join(out_dir, "motif_cooccurrence.png")
    plt.savefig(output_path, dpi=300)
    plt.close(fig)


LOGO_ALPHABET = 'ACGT'
LOGO_COLORS = {"A": '#109648', "C": '#255C99', "G": '#F7B32B', "T": '#D62839'}
LOGO_FONT = FontProperties(weight="bold")

def plot_cwms(cwms, out_dir, alphabet=LOGO_ALPHABET, colors=LOGO_COLORS, font=LOGO_FONT):
    for m, v in cwms.items():
        motif_dir = os.path.join(out_dir, m)
        os.makedirs(motif_dir, exist_ok=True)
        for cwm_type, cwm in v.items():
            output_path = os.path.join(motif_dir, f"{cwm_type}.png")

            fig, ax = plt.subplots(figsize=(10,2))

            evaluation.plot_logo(ax, cwm, alphabet, colors=colors, font_props=font)

            for name, spine in ax.spines.items():
                spine.set_visible(False)
            
            plt.savefig(output_path, dpi=100)
            plt.close(fig)




def render_image_tag(out_dir, motif, motif_type):
    """
    Renders a table cell element with an image tag or "n/a".
    
    Args:
    image_path (str): Path to the image file.
    
    Returns:
    str: A complete table cell element with the image or "n/a".
    """
    
    image_path_long = f"{out_dir}/CWMs/{motif}/{motif_type}.png"
    image_path = f"CWMs/{motif}/{motif_type}.png"
    
    if os.path.exists(image_path_long):
        return f'<img src="{image_path}" width="240">'
    else:
        return f"n/a"





def main(args):

    out_dir = args.out_dir

    os.makedirs(out_dir, exist_ok=True)

    # set up logging
    logger = logging.getLogger()
    logger.setLevel(logging.INFO)
    
    # create console handler
    ch = logging.StreamHandler()
    ch.setLevel(logging.INFO)
    
    # create formatter
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    
    # add formatter to ch
    ch.setFormatter(formatter)
    
    # add ch to logger
    logger.addHandler(ch)
    
    logger.info('loading data...')
    # load modisco objects
    modisco_obj = h5py.File(args.unified_modisco_h5)

    # load sequences in input regions and their contributions from compressed npz format
    sequences, contribs = data_io.load_regions_npz(args.regions)

    # get regions using the contribs & sequences
    if len(contribs.shape) == 3:
        regions = contribs * sequences
    elif len(contribs.shape) == 2:
        regions = contribs[:,None,:] * sequences

    # load peaks and hit calls
    half_width = regions.shape[2] // 2
    logger.info(f"half width: {half_width}")
    modisco_half_width = args.modisco_region_width // 2
    logger.info(f"modisco half width: {modisco_half_width}")

    peaks_df = data_io.load_peaks(args.peaks, None, half_width)
    hits_df = data_io.load_hits(args.hits, lazy=True)
    
    motifs_df, cwms_modisco, trim_masks, motif_names = data_io.load_modisco_motifs(
        modisco_h5_path = args.unified_modisco_h5,
        trim_threshold = 0.3,
        motif_type = "cwm",
        motifs_include = None,
        motif_name_map = None,
        motif_alphas = None,
        motif_alpha_default = args.alpha,
        include_rc = True)
    motif_names = motifs_df.filter(pl.col("motif_strand") == "+").get_column("motif_name").to_numpy()
    motif_width = cwms_modisco.shape[2]
    print(motifs_df.head())
    print(motifs_df.shape)

    # get motif names that were actually called in this dataset
    hits_df_collected = hits_df
    hits_df_collected = hits_df_collected.collect()
    
    # check which motif_names exist in the DataFrame
    motif_names_called = hits_df_collected.select(pl.col("motif_name")).unique().to_series().to_list()
    logger.info(f"number of motifs w/ >= 1 hit: {len(motif_names_called)}")

    # calculate occurrences, which is an indicator of hits per peak
    logger.info("calculating occurrences...")
    occ_df = (
            hits_df_collected
            .pivot(index="peak_id", columns="motif_name", values="count", aggregate_function="sum")
            .fill_null(0)
            .with_columns(total=pl.sum_horizontal(*motif_names_called))
            .sort(["peak_id"])
        )
    num_peaks = occ_df.height
    num_motifs = len(motif_names_called)

    logger.info("plotting hit distribution...")
    plot_hit_dist(occ_df, out_dir)

    # load modisco motif-related files
    merged_map = pl.read_csv(args.merged_motif_map, separator="\t")
    print(merged_map.shape)

    merged_anno = (
        pl.read_csv(args.merged_motif_anno, separator="\t")
        .with_columns(pl.col("pattern").alias("merged_pattern"))
        .select(["merged_pattern", "annotation"])
    )
    print(merged_anno.shape)

    map_unified_distinct, map_per_celltype = get_unified_motif_map(merged_map, merged_anno, args.celltype)
    
    # load seqlets
    logger.info("loading seqlets...")
    seqlets_df = data_io.load_modisco_seqlets(args.celltype_modisco_h5, peaks_df, half_width, modisco_half_width, lazy=True)

    # harmonize motif names: for each pattern component, find its corresponding merged pattern, and replace the motif name if that pattern exists
    seqlets_df = (
        seqlets_df
        .with_columns(
        pl.col("motif_name").map_elements(
            lambda s: map_per_celltype.filter(pl.col("component_pattern_short") == s).select(pl.col("merged_pattern")).to_series().to_list().pop(0)
                if len(map_per_celltype.filter(pl.col("component_pattern_short") == s).select(pl.col("merged_pattern")).to_series().to_list()) > 0 else "NA"
        ))
        .filter(pl.col("motif_name") != "NA")
    )

    # join hits with peaks
    hits_df_joined = (
            hits_df
            .with_columns(pl.col('peak_id').cast(pl.UInt32))
            .join(
                peaks_df.lazy(), on="peak_id", how="inner"
            )
            .select(
                chr=pl.col("chr"),
                start_untrimmed=pl.col("start_untrimmed"),
                end_untrimmed=pl.col("end_untrimmed"),
                is_revcomp=pl.col("strand") == '-',
                motif_name=pl.col("motif_name"),
                peak_region_start=pl.col("peak_region_start"),
                peak_id=pl.col("peak_id")
            )
        )

    # get unique hits
    hits_unique = hits_df_joined.unique(subset=["chr", "start_untrimmed", "motif_name", "is_revcomp"])

    # get hits overlapping modisco regions
    region_len = regions.shape[2]
    center = region_len / 2

    logger.info(f"region len: {region_len}")
    logger.info(f"center: {center}")
    logger.info(f"modisco half-width: {modisco_half_width}")

    hits_filtered = (
        hits_df_joined
        .filter(
            ((pl.col("start_untrimmed") - pl.col("peak_region_start")) >= (center - modisco_half_width)) 
            & ((pl.col("end_untrimmed") - pl.col("peak_region_start")) <= (center + modisco_half_width))
        )
        .unique(subset=["chr", "start_untrimmed", "motif_name", "is_revcomp"])
    )

    # shift hit coordinates
    motifs_discovered_from_celltype = set(seqlets_df.select(["motif_name"]).collect().to_series().to_list())
    print(motifs_discovered_from_celltype)
    print(len(motifs_discovered_from_celltype))

    celltype_modisco_obj = h5py.File(args.celltype_modisco_h5)
    hits_filtered_shifted = shift_hit_coords(seqlets_df, hits_filtered, map_per_celltype, motifs_discovered_from_celltype, modisco_obj, celltype_modisco_obj)

    # get recall and other QC stats
    recall_data, recall_df, cwms = get_recall(motif_names_called=motif_names_called,
                           seqlets_df=seqlets_df,
                           hits_unique=hits_unique,
                           hits_filtered_shifted=hits_filtered_shifted,
                           modisco_obj=modisco_obj,
                           motif_width=motif_width,
                           map_unified_distinct=map_unified_distinct,
                           regions=regions)

    out_path = os.path.join(out_dir, "seqlet_recall.tsv")
    recall_df.write_csv(out_path, separator="\t")
    
    # more plotting
    logger.info("plotting QC...")
    print(recall_df.head())
    plot_hits_per_motif(recall_df, motifs_discovered_from_celltype, out_dir)
    plot_hits_vs_restricted(recall_df, motifs_discovered_from_celltype, out_dir)
    plot_hits_vs_seqlets(recall_data, motifs_discovered_from_celltype, out_dir)
    logger.info("done plotting QC.")
    
    # get a filtered recall df which doesn't include low-quality patterns, for the report
    recall_df_filt = recall_df.filter(pl.col("annotation") != args.annotation_drop)
    
    # get most frequent motifs and plot their CWMs
    k = 50
    top_k_motifs = recall_df_filt.top_k(k, by = "num_hits_total")
    top_k_motif_names = top_k_motifs.select(pl.col("motif_name")).to_series().to_list()

    print(len(cwms))
    cwms_to_plot = {key: cwms[key] for key in cwms.keys() if key in top_k_motif_names}
    recall_df_filt2 = recall_df_filt.filter(pl.col("motif_name").is_in(top_k_motif_names))
    print(len(cwms_to_plot))
    print(recall_df_filt2.shape)

    logger.info("plotting CWMs...")
    cwm_dir = os.path.join(out_dir, "CWMs")
    # evaluation.plot_cwms(cwms_to_plot, cwm_dir, out_dir)
    plot_cwms(cwms_to_plot, cwm_dir)

    # plot hit distributions for top 10 motifs
    logger.info("plotting hit distributions...")
    plot_hit_distributions(occ_df, top_k_motif_names[0:10], out_dir)

    # calculate motif co-occurrence
    # logger.info("calculating motif occurrences/co-occurrences...")
    # calc_cooccurrence(num_peaks=num_peaks,
    #                   num_motifs=num_motifs,
    #                   occ_df=occ_df,
    #                   motif_names_called=motif_names_called,
    #                   top_k_motif_names=top_k_motif_names,
    #                   k=k,
    #                   out_dir=out_dir)

    # calculate how many hits from low quality motifs
    num_hits_drop = ( recall_df
        .filter(pl.col("annotation") == args.annotation_drop)
        .select(pl.col("num_hits_total"))
        .to_series().to_list() )
    num_hits_drop = sum(num_hits_drop)
    num_hits_total = recall_df.select(pl.col("num_hits_total")).to_series().to_list()
    num_hits_total = sum(num_hits_total)
    prop_drop = np.round(100*num_hits_drop / num_hits_total, 2)
    drop_str = f"{num_hits_drop} out of {num_hits_total} hits ({prop_drop}%) were to low-quality motifs (possibly overlapping high-quality motifs)."
    logger.info(drop_str)

    # write out the report.
    logger.info("writing report...")
    # load template then render
    with open(args.report, "r", encoding="utf-8") as f:
      template_str = f.read()
    
    template = Template(template_str)
    
    report = template.render(seqlet_recall_data=recall_df_filt.iter_rows(named=True),
                             out_dir=out_dir,
                             drop_str=drop_str,
                             motif_names=top_k_motif_names[0:10],
                             render_image_tag=render_image_tag,
                             n_motifs=len(motif_names_called))
    
    out_path = os.path.join(out_dir, "report.html")
    with open(out_path, "w") as f:
        f.write(report)

    logger.info("done.")
    

    
if __name__ == '__main__':
    args = parse_args()
    main(args)




