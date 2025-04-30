# Author: Selin Jessa
# 2024
# Purpose: helper and wrapper functions for exploring/visualizing model outputs,
# motifs, seqlets, hits, etc.

import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import logomaker


LOGO_COLORS = {"A": '#109648', "C": '#255C99', "G": '#F7B32B', "T": '#D62839'}


# def plot_logo(X, ax=None, ylab=None, title=None, xticks_fontsize=10, **kwargs):
#     """
#     Plot a sequence logo for a given sequence matrix X, which is a 2D numpy array
#     that could represent a CWM or PPM, or other similar matrix.
    
#     Args:
#         X (np.array): a 2D numpy array representing a sequence matrix (PWM, CWM)
#         ax (matplotlib.axes.Axes): the axes where the logo will be plotted
#     """
    
#     df = pd.DataFrame(X, columns=["A", "C", "G", "T"])
#     logo = logomaker.Logo(df, ax=ax, **kwargs, color_scheme=LOGO_COLORS, font_weight="bold")
#     logo.style_spines(visible=False)
#     logo.style_xticks(fmt='%d', anchor=0, fontsize=xticks_fontsize)
    
#     if ylab is not None:
#         logo.ax.set_ylabel(ylab, rotation=90, labelpad=-1)

#     if title is not None:
#         logo.ax.set_title(title, fontsize=12)
    
#     return logo




# helper function for printing a motif as a dataframe
def df(mat):
	return pd.DataFrame(mat, columns=["A", "C", "G", "T"])




LOGO_COLORS = {"A": '#109648', "C": '#255C99', "G": '#F7B32B', "T": '#D62839'}

def plot_logo(mat, title=None, type="CWM", ax=None, plot_IC=False, xticks_fontsize=10, **kwargs):
	"""
	Yet another function to plot a motif logo with some customizations.
	TODO: consolidate these.

	Args:
		mat (np.array): The motif matrix, with shape (4, L) where L is the motif length.
	"""

	if ax is None:
		plt.figure(figsize=(6, 2))
		ax = plt.subplot(111)

	# do IC-scaling if type is PWM
	if type == "PPM" and plot_IC:
		mat_df = logomaker.transform_matrix(df(mat.T), from_type='probability', to_type='information')
		type = "PWM IC (bits)"
	else:
		mat_df = pd.DataFrame(mat.T, columns=["A", "C", "G", "T"])
	
	# tangermeme.plot.plot_logo(mat, ax=ax)
	logo = logomaker.Logo(mat_df, ax=ax, **kwargs, color_scheme=LOGO_COLORS)
	logo.style_spines(visible=False)
	logo.style_xticks(fmt='%d', anchor=0, fontsize=xticks_fontsize)
	
	if title is not None:
		logo.ax.set_title(title, fontsize=12)
		 
	ax.set_ylabel(type)
	plt.tight_layout()
	
	if ax is None:
		plt.show()




def seqlet_heatmap(seqlet_ohe, pattern_name=None, width=5, height=3):
    """
    Plot a heatmap showing seqlet composition of a given pattern, by plotting
    the exact sequence of each seqlet in a pattern, colored by nucleotides. This is
    useful because it can help visualize how similar the constituent seqlets actually
    are to the consensus/average of the pattern.
    Adapted from Surag Nair:
    https://github.com/kundajelab/scATAC-reprog/blob/e93067f761eb415d772a82bddb4d714dbd732105/src/figures_factory/Fig5/high_OSK_vs_iPSC_modisco.ipynb
    """
    
    palette_nucleotides = sns.color_palette(['#109648', '#255C99', '#F7B32B', '#D62839']) 
    
    f, ax = plt.subplots(figsize=(width, height))
    
    cur = sns.heatmap(np.argmax(seqlet_ohe, -1),
            cmap = palette_nucleotides,
            cbar = 0, ax = ax)
    
    cur.set(xlabel = "position", ylabel = "seqlet")
    cur.set(xticklabels=range(0,seqlet_ohe[0].shape[0]))
    cur.set(yticklabels=[])
    # cur.set(xticklabels=[])
    # cur.set(xticks=[])
    cur.set(yticks=[])
    ax.set_title(f"{pattern_name} seqlet composition")
    
    plt.show()


def seqlet_heatmap2(seqlet_ohe, seqlet_shaps, pattern_name=None, figsize=(5,3), fade_probabilities=True, seqlet_order=None):
    """
    Plot a heatmap showing seqlet composition of a given pattern, by plotting
    the exact sequence of each seqlet in a pattern, colored by nucleotides. This is
    useful because it can help visualize how similar the constituent seqlets actually
    are to the consensus/average of the pattern. Legends and ability to fade nucleotides
    based on position importance, using plotnine.
    """
    
    from plotnine import (
        ggplot,
        aes,
        geom_tile,
        geom_text,
        scale_y_reverse,
        scale_y_discrete,
        scale_fill_brewer,
        scale_fill_manual,
        scale_alpha_continuous,
        coord_equal,
        theme,
        theme_void,
        theme_minimal,
        element_blank,
        element_rect,
        element_text,
        xlab,
        ylab,
        labs,
        ggtitle
    )

    # if an ordering of seqlets is provided (by their indices), then use those to
    # sort the OHE/shaps for visualization
    if seqlet_order is not None:
        seqlet_ohe = [seqlet_ohe[i] for i in seqlet_order]
        seqlet_shaps = [seqlet_shaps[i] for i in seqlet_order]
    
    # get a data frame containing nucleotides & contrib scores
    seqlet_ohe_df = pd.DataFrame(np.argmax(seqlet_ohe, -1))
    seqlet_ohe_df['seqlet'] = range(0, len(seqlet_ohe_df))
    # seqlet_ohe_df.head()
    
    # convert to long format
    seqlet_ohe_df_long = pd.melt(seqlet_ohe_df, id_vars = "seqlet", var_name = "position", value_name = "base")
    
    # replace positions with bases to get a nice legend
    seqlet_ohe_df_long = seqlet_ohe_df_long.replace({'base': {0: "A", 1: "C", 2: "G", 3: "T"}})
       
    # plot the heatmap, optionally adjust alpha based on loci importance
    if fade_probabilities:
        # calculate position importances within each seqlet, such that the sum of importance
        # of each seqlet is 1. i.e. how much does each position contribute to the total
        # contribution scores of the region. Credit to Salil for this formulation.
    
        # first, we only keep positive contributions
        seqlet_shaps_pos = [np.where(seq_shap < 0, 0, seq_shap) for seq_shap in seqlet_shaps]
        
        # then, drop one dimension when taking the sum:
        seqlet_imp = np.sum(seqlet_shaps_pos, axis=2, keepdims=False)
        # print(seqlet_imp.shape)
        
        seqlet_imp /= np.sum(seqlet_imp, axis=1, keepdims=True)
        # print(seqlet_imp.shape)
    
        seqlet_imp_df = pd.DataFrame(seqlet_imp.tolist())
        seqlet_imp_df['seqlet'] = range(0, len(seqlet_imp_df))
        
        # convert to long format
        seqlet_imp_df_long = pd.melt(seqlet_imp_df, id_vars = "seqlet", var_name = "position", value_name = "importance")
        
        # join dataframes containing nucleotides & importances
        seqlet_heatmap_df_long_merged = seqlet_ohe_df_long.merge(seqlet_imp_df_long, on=['seqlet', 'position'])
    
        # get the max importance to set the scale
        max_imp = max(seqlet_heatmap_df_long_merged['importance'])
        print(f"@ max importance across positions/seqlets: {max_imp}")
        
        p = (
            ggplot(seqlet_heatmap_df_long_merged, aes("position", "seqlet"))
            + geom_tile(aes(fill="factor(base)", alpha = "importance"))
            + scale_alpha_continuous(range = [0, 1], limits = [0, max_imp])
            + labs(fill = "nucleotide", alpha = "importance") + ggtitle(f"{pattern_name} seqlet composition")
        )
        
    else:
        p = (
            ggplot(seqlet_ohe_df_long, aes("position", "seqlet"))
            + geom_tile(aes(fill="factor(base)"))
            + labs(fill = "nucleotide") + ggtitle(f"{pattern_name} seqlet composition")
        )
    
    # more aesthetic adjustments
    p = (p
        + scale_fill_manual(values = {"A": '#109648', "C": '#255C99', "G": '#F7B32B', "T": '#D62839'})
        + theme( 
            axis_ticks=element_blank(),
            panel_background=element_rect(fill="white"),
            figure_size=figsize,
            legend_key_size = 3
        )
        + xlab("position") + ylab("seqlet")
    )
        
    return p





def hit_heatmap2(hit_ohe, pattern_name=None, figsize=(5,3), hit_order=None, sort_by_pos=None):
    """
    Similar to seqlet_heatmap2, but for hits. Plot a heatmap showing sequence composition
    of hits for a given pattern, by plotting the exact sequence of each hit for a pattern,
    colored by nucleotides.
    """

    from plotnine import (
        ggplot,
        aes,
        geom_tile,
        geom_text,
        geom_boxplot,
        scale_y_reverse,
        scale_y_discrete,
        scale_fill_cmap,
        scale_fill_brewer,
        scale_fill_manual,
        scale_alpha_continuous,
        scale_fill_gradient2,
        coord_equal,
        coord_cartesian,
        theme,
        theme_void,
        theme_minimal,
        element_blank,
        element_rect,
        element_text,
        xlab,
        ylab,
        labs,
        ggtitle
    )
    
    n_hits = len(hit_ohe)

    # sorting options:
    # 1. if a position to sort by is provied, use that
    # 2. if an ordering of seqlets is provided (by their indices), then use those to
    # sort the OHE/shaps for visualization
    if sort_by_pos is not None:        
        # get the nucleotide at the provided in each position, then get indices for sorting
        hit_order_pos = np.argsort([np.argmax(seq[sort_by_pos], -1) for seq in hit_ohe])
        hit_ohe = [hit_ohe[i] for i in hit_order_pos]
    
    elif hit_order is not None:
        hit_ohe = [hit_ohe[i] for i in hit_order]
    
    # get a data frame containing nucleotides & contrib scores
    hit_ohe_df = pd.DataFrame(np.argmax(hit_ohe, -1))
    hit_ohe_df['hit'] = range(0, len(hit_ohe_df))
    # seqlet_ohe_df.head()
    
    # convert to long format
    hit_ohe_df_long = pd.melt(hit_ohe_df, id_vars = "hit", var_name = "position", value_name = "base")
    
    # replace positions with bases to get a nice legend
    hit_ohe_df_long = hit_ohe_df_long.replace({'base': {0: "A", 1: "C", 2: "G", 3: "T"}})
       
    p = (
        ggplot(hit_ohe_df_long, aes("position", "hit"))
        + geom_tile(aes(fill="factor(base)"))
        + labs(fill = "nucleotide") + ggtitle(f"{pattern_name} \nhits sequence composition")
        )
    
    # more aesthetic adjustments
    p = (p
        + scale_fill_manual(values = {"A": '#109648', "C": '#255C99', "G": '#F7B32B', "T": '#D62839'})
        + theme( 
            axis_ticks=element_blank(),
            panel_background=element_rect(fill="white"),
            figure_size=figsize,
            legend_key_size = 3
        )
        + xlab("position") + ylab(f"hits (n={n_hits})")
    )
        
    return p




def normalize_row(row, range = [0, 1]):
    """
    Normalize the row to the range [0, 1] or [-1, 1].
    """

    if (range == [0, 1]):

        pseudocount = 1e-8
    
        min_val = np.min(row)
        max_val = np.max(row)
    
        return (row - min_val) / (pseudocount + max_val - min_val)

    elif (range == [-1, 1]):
        
        max_val = np.max(np.abs(row))
        return row / max_val



def trim_to_range(row, range=[-1, 1]):
    """
    Trims the values in the row to [-1, 1].
    """

    if len(row) == 0:
        return []  # Handle empty rows gracefully

    if range == [-1, 1]:
        return np.clip(row, -1, 1)

    else:
        raise ValueError("Invalid range specified. Must be [-1, 1].")


def hit_contrib_heatmap(contribution_scores, pattern_name=None, normalize=True, figsize=(10, 6), hit_order=None):
    """
    Plot a heatmap showing contribution scores at each sequence, for a set of aligned hits.
    """

    from plotnine import (
        ggplot,
        aes,
        geom_tile,
        geom_text,
        geom_boxplot,
        scale_y_reverse,
        scale_y_discrete,
        scale_fill_cmap,
        scale_fill_brewer,
        scale_fill_manual,
        scale_alpha_continuous,
        scale_fill_gradient2,
        coord_equal,
        coord_cartesian,
        theme,
        theme_void,
        theme_minimal,
        element_blank,
        element_rect,
        element_text,
        xlab,
        ylab,
        labs,
        ggtitle
    )
    
    # get order
    if hit_order is not None:
        contribution_scores = [contribution_scores[i] for i in hit_order]
        
    contribution_scores = np.array(contribution_scores)

    # normalize each row to [-1,1]
    if normalize:
        normalized_scores = np.apply_along_axis(normalize_row, 1, contribution_scores, range = [-1, 1])
    else:
        normalized_scores = contribution_scores
    
    # convert the list of arrays to a DataFrame for visualization
    score_df = pd.DataFrame(normalized_scores)
    
    n_hits = len(contribution_scores)
        
    # convert to long format
    score_df_long = pd.melt(score_df.reset_index(), id_vars = "index", var_name = "position", value_name = "score")
    score_df_long = score_df_long.rename(columns={'index': 'hit'})

    p = (
        ggplot(score_df_long, aes("position", "hit"))
        + geom_tile(aes(fill="score"))
        + labs(fill="SHAP score") + ggtitle(f"{pattern_name} \nhits contribution scores")
    )
    
    # more aesthetic adjustments
    p = (p
        + scale_fill_gradient2(low="blue", mid="white", high="red", midpoint=0)
        + theme( 
            axis_ticks=element_blank(),
            axis_text_x=element_text(size = 6),
            panel_background=element_rect(fill="white"),
            figure_size=figsize
        )
        # guides(fill=guide_colorbar(ticks.colour = NA)) +
        + xlab("position") + ylab(f"hits (n={n_hits})")
    )
    
    return p




def hit_conservation_heatmap(conservation_scores, pattern_name=None, normalize=False, trim=False, figsize=(10, 6), hit_order=None):
    """
    Plot a heatmap showing conservation scores at each sequence, for a set of aligned hits.
    """

    from plotnine import (
        ggplot,
        aes,
        geom_tile,
        geom_text,
        geom_boxplot,
        scale_y_reverse,
        scale_y_discrete,
        scale_fill_cmap,
        scale_fill_brewer,
        scale_fill_manual,
        scale_alpha_continuous,
        scale_fill_gradient2,
        coord_equal,
        coord_cartesian,
        theme,
        theme_void,
        theme_minimal,
        element_blank,
        element_rect,
        element_text,
        xlab,
        ylab,
        labs,
        ggtitle
    )
    
    # get order
    if hit_order is not None:
        conservation_scores = [conservation_scores[i] for i in hit_order]
        
    conservation_scores = np.array(conservation_scores)

    # normalize each row to [-1,1]
    if normalize:
        normalized_scores = np.apply_along_axis(normalize_row, 1, conservation_scores, range = [0, 1])
    elif trim:
        normalized_scores = np.apply_along_axis(trim_to_range, 1, conservation_scores, range = [-1, 1])
    else:
        normalized_scores = conservation_scores
    
    # convert the list of arrays to a DataFrame for visualization
    score_df = pd.DataFrame(normalized_scores)
    
    n_hits = len(conservation_scores)
        
    # convert to long format
    score_df_long = pd.melt(score_df.reset_index(), id_vars = "index", var_name = "position", value_name = "score")
    score_df_long = score_df_long.rename(columns={'index': 'hit'})

    color_limits = [np.max(score_df_long['score']), np.min(score_df_long['score'])]

    p = (
        ggplot(score_df_long, aes("position", "hit"))
        + geom_tile(aes(fill="score"))
        + labs(fill="Conservation") + ggtitle(f"{pattern_name} \nhits conservation scores")
    )
    
    # more aesthetic adjustments
    p = (p
        # + scale_fill_gradient2(low="blue", mid="white", high="red", midpoint=0, limits = color_limits)
        + scale_fill_cmap("viridis")
        + theme( 
            axis_ticks=element_blank(),
            axis_text_x=element_text(size = 6),
            panel_background=element_rect(fill="white"),
            figure_size=figsize
        )
        # guides(fill=guide_colorbar(ticks.colour = NA)) +
        + xlab("position") + ylab(f"hits (n={n_hits})")
    )
    
    return p



def hit_conservation_boxplot(conservation_scores, pattern_name=None, normalize=True,
                             figsize=(6, 3), hit_order=None):
    """
    Plot a boxplot showing the distribution of conservation scores at each position in the sequence,
    for a set of aligned hits.
    """

    from plotnine import (
        ggplot,
        aes,
        geom_tile,
        geom_text,
        geom_boxplot,
        scale_y_reverse,
        scale_y_discrete,
        scale_fill_cmap,
        scale_fill_brewer,
        scale_fill_manual,
        scale_alpha_continuous,
        scale_fill_gradient2,
        coord_equal,
        coord_cartesian,
        theme,
        theme_void,
        theme_minimal,
        element_blank,
        element_rect,
        element_text,
        xlab,
        ylab,
        labs,
        ggtitle
    )
    
    # normalize each row to [-1,1]
    if normalize:
            normalized_scores = np.apply_along_axis(normalize_row, 1, conservation_scores, range = [0, 1])
    else:
        normalized_scores = conservation_scores
    
    # convert the list of arrays to a DataFrame for visualization
    score_df = pd.DataFrame(normalized_scores)
    
    n_hits = len(conservation_scores)
        
    # convert to long format
    score_df_long = pd.melt(score_df.reset_index(), id_vars = "index", var_name = "position", value_name = "score")
    score_df_long = score_df_long.rename(columns={'index': 'hit'})

    color_limits = [np.max(score_df_long['score']), np.min(score_df_long['score'])]
    
    # get median per group
    score_df_long['score_median'] = score_df_long.groupby('position')['score'].transform('median')
    
    p = (
        ggplot(score_df_long, aes("position", "score"))
        + geom_boxplot(aes(fill="score_median"), outlier_alpha = 0.8, outlier_size = 1)
        + labs(fill="Median \nconservation")
        + ggtitle(f"{pattern_name} \nhits conservation scores")
    )
    
    # more aesthetic adjustments
    p = (p
        # + scale_fill_gradient2(low="blue", mid="white", high="red", midpoint=0, limits = color_limits)
        + scale_fill_cmap("viridis")
        + theme( 
            axis_ticks=element_blank(),
            axis_text_x=element_text(size = 6),
            panel_background=element_rect(fill="white"),
            figure_size=figsize
        )
        + coord_cartesian(ylim = [-1, 1])
        + xlab("position") + ylab(f"phyloP \nconservation score")
    )
    
    return p
