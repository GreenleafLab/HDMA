# Purpose: Merge multiple modisco patterns together, by collapsing similar patterns.
# Given a cluster key file from gimme cluster, which specifies the broad cluster,
# and the component patterns, we load the modisco h5 files for each component celltype,
# put the pattern information together, and then collapse similar patterns using
# the functionality from tfmodisco-lite. We save a TSV report that describes which
# input patterns were merged into the output/merged patterns. We save an HTML report
# that contains the info in the TSV, plus the CWM logos for the input & output/merged patterns.

# Adapted from Avanti, Ziwei, Ryan, Soumya.
# https://github.com/jmschrei/tfmodisco-lite/blob/main/examples/MergeMotifsAcrossRuns.ipynb

from modiscolite.aggregator import SimilarPatternsCollapser
from modiscolite.core import TrackSet, Seqlet, SeqletSet
from modiscolite.report import _plot_weights,path_to_image_html
import h5py as h5
import numpy as np
import logomaker
import os
import matplotlib.pyplot as plt
import pandas as pd
from datetime import datetime
import argparse

# need full paths!
def parse_args():
    parser = argparse.ArgumentParser(description="Merge modisco patterns.", formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("--out-dir", type=str, required=True, help="Output directory, specifying where to write collapsed .h5 files, and summary report.")
    parser.add_argument("--model-head", type=str, required=True,
                       help='''
                       Either 'counts' or 'profile'. Indicates whether contribs were derived with respect
                       to counts or profile model head. Used for finindg correct modisco/contribs h5 files.
                       ''')
    parser.add_argument("--cluster-key", type=str, required=True,
                        help='''
                        cluster_key.txt containing output of gimme cluster.
                        TSV file without header expected to contain at least four columns:
                        1) gimme cluster id (denoting a group of patterns to merge),
                        2) pattern class (either 'pos_patterns' or 'neg_patterns'),
                        3) component patterns as a comma separated list of all the component patterns,
                        and 4) batch (so that this script can operate on a subset of clusters
                        to aid parallelization).
                        Each component pattern in the 3rd column should be in the form
                        [celltype]__[pattern_class].[pattern_id], e.g. 'Brain_c1__pos_patterns.pattern_1'.
                        Note the double underscores. Every [celltype] is expected to have a modisco h5 file
                        and contribs h5 file.
                        ''')
    parser.add_argument("--modisco-dir", type=str, required=True,
                        help='''
                        Directory holding one directory per [celltype], each of which contains
                        a modisco h5 file, named as [model_head]_modisco_output.h5.
                        ''')
    parser.add_argument("--contribs-dir", type=str, required=True,
                        help=''' 
                        Directory holding one directory per [celltype], each which contains
                        a contribs h5 file, containing contribution scores averaged
                        across folds, named as average_shaps.[model_head].h5.
                        ''')
    parser.add_argument("--batch", type=int, required=True, help="Batch number to process.")
    parser.add_argument("--debug", type=bool, required=False, default=False)

    args = parser.parse_args()
    print(args)
    
    # check that the model head is valid
    assert args.model_head in ["counts", "profile"]

    # check that the cluster key file exists
    assert os.path.exists(args.cluster_key), f"Cluster key file {args.cluster_key} not found."

    # check that the modisco directory exists
    assert os.path.exists(args.modisco_dir), f"Modisco directory {args.modisco_dir} not found."

    # check that the contribs directory exists
    assert os.path.exists(args.contribs_dir), f"Contribs directory {args.contribs_dir} not found."

    # check that the output directory exists
    assert os.path.exists(args.out_dir), f"Output directory {args.out_dir} not found."

    return args


def main(args):
    
    # SET UP ARGS -----------------------------------------------------------------
    model_head=args.model_head
    modisco_dir=args.modisco_dir
    contribs_dir=args.contribs_dir
    out_dir=args.out_dir
    
    
    # SET DEFAULTS ----------------------------------------------------------------
    
    # tfmodisco-lite defaults (as of version 2.2.0 from https://github.com/jmschrei/tfmodisco-lite/blob/main/examples/MergeMotifsAcrossRuns.ipynb)
    min_overlap  = 0.7
    prob_and_pertrack_sim_merge_thresholds = [(0.8,0.8), (0.5, 0.85), (0.2, 0.9)]
    prob_and_pertrack_sim_dealbreaker_thresholds = [(0.4, 0.75), (0.2,0.8), (0.1, 0.85), (0.0,0.9)]
    min_frac = 0.2 # also called frac_support_to_trim_to
    min_num = 30 # also called min_num_to_trim_to
    flank_to_add = 5 # also called initial_flank_to_add
    window_size = 20 # also called trim_to_window_size
    max_seqlets_subsample = 300 # also called merging_max_seqlets_subsample
    # track_signs = 1 for positive, -1 for negative
    input_window = 400 # corresponds to default for tfmodisco-lite args.window
    def calculate_window_offsets(center: int, window_size: int) -> tuple:
        return (center - window_size // 2, center + window_size // 2)
    
    
    
    # GET PATTERNS -----------------------------------------------------------------
    
    # loop through all gimme-clusters and extract the component patterns
    print("@ reading cluster key file...")
    pattern_clusters_to_merge = []
    
    with open(args.cluster_key) as clusters:
        for line in clusters:
            line = line.strip()
            fields = line.split("\t")
            cluster_id = fields[0] # the "id" of the gimme cluster
            pattern_class = fields[1] # either pos_patterns or neg_patterns
            component_patterns = fields[2].split(",") # list of component patterns in the form [celltype]__[pattern_class].[pattern_id]
            batch = fields[3] # batch number
            components = []
            
            # for each component pattern, extract the celltype and pattern name
            for pattern in component_patterns:
                pattern_components = pattern.split("__")
                pattern_celltype = pattern_components[0]
                # split off the positive/negative prefix, because patterns names don't have those
                # prefixes in the modisco h5 files
                pattern_name = pattern_components[1].split(".")[1]
                components.append({ "celltype": pattern_celltype,
                                    "name": pattern_name
                })
            pattern_clusters_to_merge.append({
                "id": cluster_id,
                "pattern_class": pattern_class,
                "components": components,
                "batch": batch
            })
    
    if args.debug:
        # just get one
        print("@ debugging with first cluster.")
        pattern_clusters_to_merge = [pattern_clusters_to_merge[0]]
    
    elif args.batch is not None:
        # filter to only the batch we want to process
        print(f"@ filtering clusters to batch {args.batch}")
        pattern_clusters_to_merge = [cluster for cluster in pattern_clusters_to_merge
                                     if cluster["batch"] == str(args.batch)]
                                     
   # [print(cluster["batch"]) for cluster in pattern_clusters_to_merge]
                                     
    # print(pattern_clusters_to_merge)
    print("@ processing " + str(len(pattern_clusters_to_merge)) + " clusters.")
    
    # get unique cell types and their associated patterns that belong to this pattern cluster
    for pattern_cluster in pattern_clusters_to_merge:
    
        print("@ starting cluster " + pattern_cluster['id'])
        print(pattern_cluster)

        cluster_out_dir = out_dir + "/" + pattern_cluster['id']
        os.makedirs(cluster_out_dir, exist_ok=True)

        pattern_class = pattern_cluster["pattern_class"]
        print("\t@ pattern class: " + pattern_class)
        
        if len(pattern_cluster["components"]) == 1:
            print("\t@ found cluster with only 1 component pattern.")
        
        # set up output paths
        merged_modisco_h5_path = f"{cluster_out_dir}/merged_modisco.h5"
        merged_modisco_tsv = f"{cluster_out_dir}/merged_modisco.tsv"
        merged_modisco_html = f"{cluster_out_dir}/merged_modisco.html"
    
        # merge lists
        patterns_to_merge = []
          
        # get unique cell types and their associated patterns that belong to this cluster
        component_celltypes = set([component["celltype"] for component in pattern_cluster["components"]])
          
        # get a dict of celltype:patterns
        for component_celltype in component_celltypes:
            component_patterns = [component["name"] for component in pattern_cluster["components"]
                                if component["celltype"] == component_celltype]
            patterns_to_merge.append((component_celltype, component_patterns))
    
        print("\t@ merging " + str(len(pattern_cluster["components"])) + " patterns...")
    
        # set up empty lists
        union_onehot = []
        union_hypscores = []
        union_projscores = []
        union_patterncoords = []	
    
        # get an index of each pattern, so we can match the unioned patterns to the original patterns
        all_pattern_idxs = []
        all_pattern_idxs = [{celltype: pattern_name} for celltype, patterns in patterns_to_merge for pattern_name in patterns]

        # LOAD MODISCO H5 FILES --------------------------------------------------------
        # for each cluster, load the modisco h5 files and extract the pattern information
        print("\t@ loading modisco h5 files...")
        # load in the modisco h5s for each cell type, and extract info 
        print(f"\t@ extracting data per celltype pattern [{datetime.now().strftime('%H:%M:%S')}]")
        exampleidx_offset = 0 # incremented after each modisco results file
    
        for celltype, patterns in patterns_to_merge:
            print("\t\t@ on celltype", celltype, "with patterns:", patterns)
    
            shaps_file = f"{contribs_dir}/{celltype}/average_shaps.{model_head}.h5"
            results_file = f"{modisco_dir}/{celltype}/{model_head}_modisco_output.h5"
    
            with h5.File(shaps_file) as shaps_fh, h5.File(results_file) as results_fh:
                onehot = shaps_fh["raw"]["seq"]
                hypscores = shaps_fh["shap"]["seq"]
                projscores = shaps_fh["projected_shap"]["seq"]
                
                allpattern_exampleidxs = []
                # first, iterate through the patterns and get all the example indices
                # (Note: "example_idx" refers to the index of the sequence that contained
                #  the seqlet)
                for pattern_name in patterns:
                    seqlets_grp = results_fh[pattern_class][pattern_name]['seqlets']
                    allpattern_exampleidxs.extend(np.array(seqlets_grp['example_idx']))
                
                # figure out the subset of indices that actually have seqlets, sort it.
                surviving_indices = sorted(list(set(allpattern_exampleidxs)))
                print("\t\t" + str(len(surviving_indices))+" indices had seqlets out of "
                    +str(len(onehot)))
                # add the scores for the subset of sequences that have scores to the
                # 'union' list.
                center = onehot.shape[2] // 2
                start, end = calculate_window_offsets(center, input_window)
                for idx in surviving_indices:
                    union_onehot.append(onehot[idx,:,start:end])
                    union_hypscores.append(hypscores[idx,:,start:end])
                    union_projscores.append(projscores[idx,:,start:end])
                
                # create an index remapping based on the subset of surviving indices
                # (we will add exampleidx_offset later)
                idx_remapping = dict(zip(surviving_indices,
                                    np.arange(len(surviving_indices))))
                
                # Now iterate through the patterns again and prep the seqlet coordinates,
                # remapping the example indices as needed.
                # We also add in exampleidx_offset to account for all the previous seqeuences
                # that have already been added to the 'union' lists
                for pattern_name in patterns:
                    seqlets_grp = results_fh[pattern_class][pattern_name]['seqlets']
                    pattern_exampleidxs = np.array(seqlets_grp['example_idx'])
                    # remap the example idxs
                    pattern_remapped_exampleidxs = np.array([
                        (exampleidx_offset+idx_remapping[idx]) for idx in pattern_exampleidxs])
                    pattern_start = np.array(seqlets_grp['start'])
                    pattern_end = np.array(seqlets_grp['end'])
                    pattern_isrevcomp = np.array(seqlets_grp['is_revcomp'])
                    union_patterncoords.append((pattern_remapped_exampleidxs,
                                                pattern_start, pattern_end, pattern_isrevcomp))
                # increment exampleidx_offset
                exampleidx_offset = (exampleidx_offset + len(surviving_indices))
        print(f"\t@ finished data per celltype/pattern [{datetime.now().strftime('%H:%M:%S')}]")
    
        # reshape the combined info and build track set
        union_onehot = np.transpose(np.array(union_onehot), axes=(0,2,1))
        union_hypscores = np.transpose(np.array(union_hypscores), axes=(0,2,1))
        union_projscores = np.transpose(np.array(union_projscores), axes=(0,2,1))
        track_set = TrackSet(one_hot=union_onehot,
                            contrib_scores=union_projscores,
                            hypothetical_contribs=union_hypscores)
    
        print("\t@ combining pattern information...")
        #Create pattern objects using the new track_set and modified coordinates
        print("\t@ creating list of all seqlets.")
        all_patterns = []
        for (example_idxs, starts, ends, isrevcomps) in union_patterncoords:
            #tfmlite reuses the same object for representing seqlet
            # coordinates as well as seqlets
            seqlet_coords = [Seqlet(example_idx, start, end, isrevcomp) for
                        (example_idx, start, end, isrevcomp) in zip(
                        example_idxs, starts, ends, isrevcomps)]
            seqlets = track_set.create_seqlets(seqlet_coords)
            pattern = SeqletSet(seqlets) #SeqletSet in tfm lite = AggregatedSeqlet in tfm
            all_patterns.append(pattern)
            print("\t\t@ adding following num of seqlets:", str(len(seqlets)))
    
        # nucleotide background frequency
        bg_freq = np.mean(union_onehot, axis=(0, 1))
        # print(bg_freq)
    
        # MERGE PATTERNS -----------------------------------------------------------
        # merge and collapse patterns
        print(f"\t@ merging patterns [{datetime.now().strftime('%H:%M:%S')}]")
        merged_patterns, pattern_merge_hierarchy = SimilarPatternsCollapser(
            patterns=all_patterns,
            track_set=track_set,
            min_overlap=min_overlap,
            prob_and_pertrack_sim_merge_thresholds=prob_and_pertrack_sim_merge_thresholds,
            prob_and_pertrack_sim_dealbreaker_thresholds=prob_and_pertrack_sim_dealbreaker_thresholds,
            min_frac=min_frac,
            min_num=min_num,
            flank_to_add=flank_to_add,
            window_size=window_size,
            bg_freq=bg_freq,
            max_seqlets_subsample=max_seqlets_subsample)
        print(f"\t@ finished merging [{datetime.now().strftime('%H:%M:%S')}]")
    
        print("\t@ found", str(len(merged_patterns)), "merged patterns.")
        merged_pattern_names = []
        for i in range(len(merged_patterns)):
    
            merged_pattern_name = f"{pattern_cluster['id']}__merged_pattern_{i}"
            merged_pattern_names.append(merged_pattern_name)
            print(f"\t\t{merged_pattern_name}, supported by {len(merged_patterns[i].seqlets)} seqlets")
    
        # SAVE MERGED PATTERNS -----------------------------------------------------
        # save to disk
    
        with h5.File(merged_modisco_h5_path, "w") as output_h5:
            
            output_h5_group = output_h5.create_group(pattern_class)
            
            for i in range(len(merged_patterns)):
    
                merged_pattern = merged_patterns[i]
                output_h5_group.create_dataset(f"{pattern_cluster['id']}__merged_pattern_{i}/contrib_scores", data=merged_pattern.contrib_scores)
                output_h5_group.create_dataset(f"{pattern_cluster['id']}__merged_pattern_{i}/sequence", data=merged_pattern.sequence)
                output_h5_group.create_dataset(f"{pattern_cluster['id']}__merged_pattern_{i}/hypothetical_contribs", data=merged_pattern.hypothetical_contribs)
            
        print("\t@ done writing to disk.")
    
        # SAVE SUMMARY REPORT ------------------------------------------------------
        # Write a report that summarizes which input patterns were merged into the
        # final patterns after collapsing. This leverages the pattern_merge_hierarchy which stores
        # a tree of which patterns were merged. The original/input patterns are the child nodes 
        # in the hierarchy.

    	# make output directory
        os.makedirs(cluster_out_dir + "/logos", exist_ok=True)

        # decide whether to make HTML report as well as TSV (between 2-100 patterns being merged)
        if len(pattern_cluster["components"]) > 1 and len(pattern_cluster["components"]) < 50:
            print("\t@ making HTML and TSV reports.")
            make_logo = True
        
        else:
            print("\t@ making TSV report.")
            make_logo = False

        # make the merged logos for all pattern clusters
        merged_logo_paths = []
        # first, make a logo for each merged pattern
        for i in range(len(merged_patterns)):
        
            merged_pattern = merged_patterns[i]
            logo_path = f"logos/{pattern_cluster['id']}__merged_pattern_{i}.png"
            _plot_weights(merged_pattern.contrib_scores, cluster_out_dir + "/" + logo_path)
            merged_logo_paths.append("./" + logo_path)

        def find_pattern_matches(node, merged_pattern_name, all_patterns_idxs, merged_logo_path, make_logo = False):
    
            matched_indices = np.nonzero([np.allclose(node.pattern.contrib_scores, x.contrib_scores) for x in all_patterns])[0]
            assert len(matched_indices)==1
            # print(f"\tmatches index {matched_indices[0]} in the original patterns: {all_pattern_idxs[matched_indices[0]]}")
    
            source_celltype = next(iter(all_pattern_idxs[matched_indices[0]].keys()))
            source_pattern = next(iter(all_pattern_idxs[matched_indices[0]].values()))
            
            if make_logo:

                # produce forward and reverse logos for the input CWM
                logo_path_fwd = f"logos/{source_celltype}__{pattern_class}.{source_pattern}.fwd.png"
                _plot_weights(node.pattern.contrib_scores, cluster_out_dir + "/" + logo_path_fwd)

                cwm_rev = node.pattern.contrib_scores[::-1, ::-1]
                logo_path_rev = f"logos/{source_celltype}__{pattern_class}.{source_pattern}.rev.png"
                _plot_weights(cwm_rev, cluster_out_dir + "/" + logo_path_rev)

                logo_path_fwd = "./" + logo_path_fwd
                logo_path_rev = "./" + logo_path_rev
        
            else:
                logo_path_fwd = None
                logo_path_rev = None
    
            return [merged_pattern_name, source_celltype, source_pattern, merged_logo_path, logo_path_fwd, logo_path_rev]
    
        def pattern_match_dfs(node, merged_pattern_name, merged_logo_path, make_logo = False):
            """
            Performs a depth-first search on the PatternMergeHierarchy structure.
            
            Args:
            node: The starting node for the search.
            """
    
            results = []
            
            # print("\tnum seqlets", len(node.pattern.seqlets), "futher child nodes:", child_node.child_nodes)
            
            # check if the node has child nodes
            if len(node.child_nodes) > 0:
                # recursively visit each child node
                for child in node.child_nodes:
                    child_results = pattern_match_dfs(child, merged_pattern_name, merged_logo_path, make_logo = make_logo)
                    results.extend(child_results)
            else:
                # we're at a leaf node (empty child_nodes list)
                outs = find_pattern_matches(node, merged_pattern_name, all_pattern_idxs, merged_logo_path, make_logo = make_logo)
                results.append({'merged_pattern': outs[0],
                                'merged_logo': outs[3],
                                'input_logo_fwd': outs[4],
                                'input_logo_rev': outs[5],
                                'component_celltype': outs[1],
                                'pattern_class': pattern_class,
                                'pattern': outs[2],
                                'n_seqlets': len(node.pattern.seqlets)})
    
            return results if results else None
    
        # combine results from all root patterns
        all_results = pd.concat(
            [pd.DataFrame(pattern_match_dfs(root, merged_pattern_names[i], merged_logo_paths[i], make_logo = make_logo))
             for i, root in enumerate(pattern_merge_hierarchy.root_nodes)],
            ignore_index=True)

        # save HTML
        if make_logo:
            all_results.to_html(open(merged_modisco_html, 'w'),
            		escape=False, formatters=dict(merged_logo=path_to_image_html, input_logo_fwd=path_to_image_html, input_logo_rev=path_to_image_html), 
            		index=False)

        # save TSV
        all_results.to_csv(merged_modisco_tsv, sep="\t", index=False, columns = ["merged_pattern", "component_celltype", "pattern_class", "pattern", "n_seqlets"])

        print("\t@ wrote summary report(s).")
        print("@ done with cluster ", pattern_cluster['id'])
    
    # create a file indicating we're done.
    def touch(fname): open(fname, 'w').close()
    
    if not args.debug:
        touch(f"{out_dir}/.{args.batch}.done")


if __name__=="__main__":
    args = parse_args()
    main(args)
