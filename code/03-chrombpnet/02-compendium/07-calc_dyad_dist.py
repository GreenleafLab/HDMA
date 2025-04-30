# Purpose: given nucleosome dyad calls from NucleoATAC, calculate the distribution
# of pairwise distances between dyads.

import pandas as pd
import numpy as np
from scipy.spatial import cKDTree
import argparse


def read_bed(file_path):
    """
    Read BED file w/ nucleosome positions.
    """
  
    # read the first three columns (chrom, start, end) from the BED file, ignore the rest
    bed_columns = ['chrom', 'start', 'end']
    bed_df = pd.read_csv(file_path, sep='\t', usecols=[0, 1, 2], header=None, names=bed_columns)
    return bed_df


def pairwise_distances_per_region(bed_df, max_dist=250, k_nearest=100):
    """
    Calculate pairwise distances for each region.
    """
    
    # since regions are length 1, we use the 'start' position as the coordinate
    bed_df['position'] = bed_df['start']
    
    # separate the BED data by chromosome to avoid cross-chromosomal calculations
    all_distances = []

    for chrom, chrom_df in bed_df.groupby('chrom'):
      
        print(chrom)
      
        # create KDTree for the regions on this chromosome
        coordinates = chrom_df['position'].values.reshape(-1, 1)
        tree = cKDTree(coordinates)

        # iterate over each region to find its nearest neighbors
        for i in range(len(chrom_df)):
            
            # query the nearest k_nearest regions
            distances, indices = tree.query(coordinates[i], k=k_nearest+1)  # k+1 to exclude itself
            
            # remove distance to itself (which is 0)
            distances = distances[distances > 0]
            
            # filter distances below max_dist
            distances = distances[distances < max_dist]
            
            if len(distances) > 0:
                all_distances.extend(distances)
    
    return all_distances


def bin_distances(distances, bin_size=10, num_bins=25):
    """
    Bin distances into 25 bins (10, 20, 30, etc.)
    """
  
    max_bin = bin_size * num_bins
    bins = np.arange(0, max_bin + bin_size, bin_size)
    binned_distances, bin_edges = np.histogram(distances, bins=bins)

    # bin labels as the upper edge of each bin (10, 20, 30, etc.)
    bin_labels = [str(int(bin_edges[i])) for i in range(1, len(bin_edges))]
    return binned_distances, bin_labels


def main():
    parser = argparse.ArgumentParser(description='Calculate pairwise distances between nucleosome dyad position calls.')
    parser.add_argument('bed', type=str, help='Path to the BED file containing genomic regions')
    parser.add_argument('out_prefix', type=str, help='Output prefix for binned distances')
    parser.add_argument('--max_dist', type=int, default=250, help='Maximum distance threshold for pairwise distances (default: 250)')
    parser.add_argument('--k_nearest', type=int, default=100, help='Number of nearest neighbors to consider (default: 100)')
    parser.add_argument('--out', type=str, help='Output file to save binned distances')
    args = parser.parse_args()

    # load BED file
    bed_df = read_bed(args.bed)
    print(bed_df.head())

    # calculate distances
    distances = pairwise_distances_per_region(bed_df, max_dist=args.max_dist, k_nearest=args.k_nearest)

    # bin distances into 25 bins (10, 20, 30, etc.)
    binned_distances, bin_labels = bin_distances(distances)

    # make a df
    bin_df = pd.DataFrame({
        'Bin': bin_labels,
        'Count': binned_distances
    })

    # save out
    bin_df.to_csv(args.out_prefix + ".dyad_binned_distances.tsv", index=False, sep = "\t")
    # print(f"Total number of distances: {len(distances)}")
    print(bin_df)
        
if __name__ == '__main__':
    main()
