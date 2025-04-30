# Purpose: after running FiNeMo-GPU for hit calling, generate a filtered and reconciled
# set, which does the following:
# - combine hits from hitcalling with possibly two sets of parameters
# - filter out hits to low-quality motifs
# - drop hits associated with patterns that have a poor correlation between hit-CWM
#    and MoDISco-CWM
# - for the remaining hits, if multiple hits to motifs with the same annotation
#   overlap by more than X bp, then keep only the hit with the highest hit score
# Finally, input motif names are replaced with annotated motif names, and output
# to a bed file.
#
# Optionally, this script can be run with --label-only to simply relabel motifs
# in the BED fole with a provided annotation so they can be viewed in a genome browser.
# No filtering will be performed.

import sys
import os
import warnings
import pandas as pd
import pyranges as pr
import numpy as np
import argparse
import logging


def parse_args():
    parser = argparse.ArgumentParser()
    
    # general arguments
    parser.add_argument("--out-dir", type=str, default=None, help="Output directory for results.")
    parser.add_argument("--merged-motif-anno", type=str, default=None,
                      help="Path to the merged motif annotation file (tsv). Expects columns 'pattern', 'annotation', 'idx_uniq'.")
    parser.add_argument("--hits", type=str, default=None,
                      help="The `hits_unique.tsv` output file from `finemo call-hits`.")
    parser.add_argument("--hits2", type=str, default=None,
                      help="Optionally, another `hits_unique.tsv` output file from a different run of `finemo call-hits`.")
    parser.add_argument("--recall", type=str, default=None,
                      help="The `seqlet_recall.tsv` output file from the reporting step.")
    parser.add_argument("--recall2", type=str, default=None,
                      help="Optionally, another `seqlet_recall.tsv` output file reporting on a different run.")
    parser.add_argument("--overlap-distance", type=int, default=3,
                        help="If multiple hits to motifs with the same annotation overlap by more than this distance, keep only the hit with the highest hit score.")
    parser.add_argument("--cor-threshold", type=float, default=0.9,
                       help="Threshold used to filter hits; *all* hits assigned to patterns where the hit-CWM to MoDISco-CWM correlation is greater than this value will be retained.")
    parser.add_argument("--reconcile-by", type=str, default="hit_correlation",
                        help="For reconciling overlapping hits with the same annotation, the column to use to choose the one to retain.")
    parser.add_argument("--annotation-drop", type=str, default="low-quality",
                        help="Hits associated with this annotation will be filtered out.")
    parser.add_argument("--label-only", action="store_true", default=False,
                        help="If True, only label hits using the provided annotation and do not filter or reconcile.")
    args = parser.parse_args()
    print(args)

    # check that the files exit
    assert os.path.exists(args.merged_motif_anno), f"File {args.merged_motif_anno} not found."
    assert os.path.exists(args.hits), f"File {args.hits} not found."

    # if the second hits file is provided, check that it exists
    if args.hits2 is not None:
        assert os.path.exists(args.hits2), f"File {args.hits2} not found."
        assert args.recall2 is not None, "Recall file for hits2 is required."
        
    if args.recall2 is not None:
       assert os.path.exists(args.recall2), f"File {args.recall2} not found."

    # check if the reconcile-by column is among valid options
    assert args.reconcile_by in ["hit_correlation", "hit_importance", "hit_coefficient"], "Invalid column for reconciling overlapping hits with the same annotaiton."

    return args




def main(args):

    out_dir = args.out_dir

    if args.label_only:
        out_path_bed = os.path.join(out_dir, "hits.labeled.bed")
    else:
        out_path_tsv = os.path.join(out_dir, "hits_unique.reconciled.tsv")
        out_path_bed = os.path.join(out_dir, "hits.reconciled.bed")
        out_path_hits_per_motif = os.path.join(out_dir, "hits_per_motif.tsv")
        out_path_hits_per_peak = os.path.join(out_dir, "hits_per_peak.tsv")
        out_path_peaks_top_hits = os.path.join(out_dir, "peaks_with_top_10_n_hits.tsv")

    half_width=500

    # set up logging
    logger = logging.getLogger()
    logger.setLevel(logging.INFO)
    ch = logging.StreamHandler()
    ch.setLevel(logging.INFO)
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    ch.setFormatter(formatter)
    logger.addHandler(ch)

    # start code
    # NOTE: modify the handling of motif names here to match a different structure as needed.
    logger.info('getting motif annotation...')
    merged_anno = pd.read_csv(args.merged_motif_anno, sep = "\t")
    merged_anno = merged_anno.rename(columns={"pattern": "merged_pattern"})
    merged_anno["pattern_class"] = merged_anno["pattern_class"].map(lambda s: "pos_patterns" if s == "pos" else "neg_patterns")
    merged_anno["merged_pattern"] = merged_anno["pattern_class"] + "." + merged_anno["merged_pattern"]
    merged_anno = merged_anno[["idx_uniq", "merged_pattern", "annotation", "pattern_class"]]

    print(merged_anno.head())

    logger.info('loading hits...')
    hits_df = pd.read_csv(args.hits, sep="\t")
    hits_df['source'] = '1'

    if args.label_only:

        logger.info('skipping filtering, labeling motifs only...')
        hits_filt = hits_df

        # merge with labels
        hits_filt_lab = hits_filt.merge(
            merged_anno.rename(columns={"merged_pattern": "motif_name"}),
            on="motif_name",
            how="left"
        )
        
        hits_filt_lab = hits_filt_lab.rename(columns={"motif_name": "motif_name_unlabeled"})
        hits_filt_lab['motif_name'] = hits_filt_lab['idx_uniq'].map(str) + '|' + hits_filt_lab['annotation']
        hits_final = hits_filt_lab.drop(columns=['annotation'])

    else:
    
        # filter out patterns with poor correlation between hit-CWM and MoDISco-CWM
        logger.info('filtering out patterns with poor CWM correlations...')
        # read in recall stats from reporting
        recall_df = pd.read_csv(args.recall, sep="\t")

        # filter out patterns where the hit CWM correlates poorly with the MoDISco CWM
        patterns_keep = recall_df[recall_df["cwm_correlation"] > args.cor_threshold]["motif_name"].tolist()
        print(len(patterns_keep))

        hits_filt = hits_df[hits_df["motif_name"].isin(patterns_keep)]
        hits_filt.shape
        hits_filt.head(20)

        # if a second hits file is provided, load it and combine with the first set of hits
        if args.hits2 is not None:

            logger.info('loading and filtering hits2...')
            hits2_df = pd.read_csv(args.hits2, sep="\t")
            hits2_df['source'] = '2'
            
            # read in recall stats from reporting
            recall2_df = pd.read_csv(args.recall2, sep="\t")
        
            # filter out patterns where the hit CWM correlates poorly with the MoDISco CWM
            patterns_keep2 = recall2_df[recall2_df["cwm_correlation"] > args.cor_threshold]["motif_name"].tolist()
            print(len(patterns_keep2))
        
            hits2_filt = hits2_df[hits2_df["motif_name"].isin(patterns_keep2)]
            hits2_filt.shape
            hits2_filt.head(20)
            
            # proceed with the concatenated, filtered hits
            hits_filt = pd.concat([hits_filt, hits2_filt])

        print(hits_filt.shape)
        print(hits_filt.head(20))

    
        # merge with labels
        hits_filt_lab = hits_filt.merge(
            merged_anno.rename(columns={"merged_pattern": "motif_name"}),
            on="motif_name",
            how="left"
        )

        hits_filt_lab = hits_filt_lab.rename(columns={"motif_name": "motif_name_unlabeled",
                                                    "chr": "Chromosome",
                                                    "start": "Start",
                                                    "end": "End"})

        # drop hits to low quality motifs
        logger.info('filtering out low-quality patterns...')
        hits_filt2_lab = hits_filt_lab[hits_filt_lab['annotation'] != args.annotation_drop]
        print(hits_filt2_lab.head(20))
        hits_filt2_lab.shape

        logger.info('find overlaps and reconcile hits...')
        # convert to PyRanges object and cluster overlapping intervals, ignoring strand
        gr = pr.PyRanges(hits_filt2_lab).cluster(slack = -args.overlap_distance, strand=None)
        print(gr.head(20))

        hits_filt2_lab_clust = gr.df
        hits_filt2_lab_clust.head()

        # within each cluster, for hits with the same annotation,
        # pick the one with the highest correlation. Here, I use correlation instead
        # of the hit coefficient b/c more directly comparable interpretation
        # for different hitcalling runs than hit coefficient.
        hits_filt2_lab_clust_rec = hits_filt2_lab_clust.loc[hits_filt2_lab_clust.groupby(['Cluster', 'annotation'])[args.reconcile_by].idxmax()]
        hits_filt2_lab_clust_rec.head(20)

        logger.info('writing reconciled hits to out...')
        # create unique 'motif_name' by concatenating 'idx_uniq', a unique index
        # per motif, and 'annotation' with '|'
        hits_filt2_lab_clust_rec['motif_name'] = hits_filt2_lab_clust_rec['idx_uniq'].map(str) + '|' + hits_filt2_lab_clust_rec['annotation']
        hits_final = hits_filt2_lab_clust_rec.drop(columns=['annotation'])

        # rename columns and reorder columns for BED format
        hits_final = hits_final.rename(columns={"Chromosome": "chr",
                                            "Start": "start",
                                            "End": "end"})
        
    # reorder columns
    hits_final = hits_final[
        ["chr", "start", "end", "start_untrimmed", "end_untrimmed", "motif_name",
        "source", "hit_coefficient", "hit_correlation", "hit_importance", "strand",
        "peak_name", "peak_id", "motif_name_unlabeled", "pattern_class"]
    ]

    # need to write bed without start trimmed values so the format is correct
    if not args.label_only:
      
        hits_final.to_csv(out_path_tsv, sep = "\t", index=False, header=True)
        
        # write out the total hits per motif
        hits_final_grouped = hits_final.groupby(['motif_name', 'pattern_class']).size().reset_index(name='count')
        hits_final_grouped.sort_values(by='count', ascending=False).to_csv(out_path_hits_per_motif, sep = "\t", index=False, header=True)
        
        # write out stats about the number of hits per peak, grouped by positive or negative patterns
        hits_per_peak = hits_final.groupby(['peak_id', 'pattern_class']).size().reset_index(name='count')
    
        # function to get stats
        def get_stats(pattern_class=None):
          
          if pattern_class:
              filtered_hits = hits_per_peak[hits_per_peak["pattern_class"] == pattern_class]['count']
              suffix = pattern_class
          else:
              suffix = "all"
              filtered_hits = hits_per_peak['count']
      
          return {
              f'min_{suffix}': [filtered_hits.min()],
              f'median_{suffix}': [filtered_hits.median()],
              f'mean_{suffix}': [filtered_hits.mean()],
              f'max_{suffix}': [filtered_hits.max()]
          }

        # get stats for both pos_patterns and neg_patterns
        hits_per_peak_stats = {
          **get_stats('pos_patterns'),
          **get_stats('neg_patterns'),
          **get_stats()}

        pd.DataFrame.from_dict(hits_per_peak_stats).to_csv(out_path_hits_per_peak, sep = "\t", index=False, header=True)
        
        # for troubleshooting, save the top 10 peaks with highest number of hits
        peaks_top_hits = hits_per_peak.sort_values(by='count', ascending=False).head(10)
        peaks_top_hits.to_csv(out_path_peaks_top_hits, sep = "\t", index=False, header=True)
        
    # multiply hit score by 1000 to get a score for viewing in browser
    hits_final['hit_correlation'] = hits_final['hit_correlation'] * 1000
    hits_bed = hits_final[["chr", "start", "end", "motif_name", "hit_correlation", "strand", "pattern_class"]]

    print(hits_bed.head(20))

    hits_bed.to_csv(out_path_bed, sep = "\t", index=False, header=False)
    logger.info("done.")


    
if __name__ == '__main__':
    args = parse_args()
    main(args)

