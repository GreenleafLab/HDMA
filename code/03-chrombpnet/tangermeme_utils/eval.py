# Author: Selin Jessa
# 2024

import torch
import numpy as np

def mean_across_sequences(y):
    """
    Compute the mean across sequences (dimension 1) for each fold.

    Parameters
    ----------
    tensor: torch.Tensor
        Input tensor where the first two dimensions are expected to be folds and sequences,
        shape (5, sequences, ...).

    Returns
    -------
    mean_seq: torch.Tensor
        Tensor of means across sequences, shape (folds, positions, ...).
    """
    y_across_seqs = torch.mean(y, dim=1)
    return y_across_seqs






def summarize_across_folds(tensor):
    """
    Compute the mean and standard deviation across folds (dimension 0).

    Parameters
    ----------
    tensor: torch.Tensor
        Tensor where the first dimension is expected to correspond to (folds, ...).

    Returns
    -------
    mean_folds: torch.Tensor
        Tensor of means across folds, shape (positions).

    sd_folds: torch.Tensor
        Tensor of standard deviations across folds, shape (positions).
    """
    
    mean_folds = torch.mean(tensor, dim=0)
    sd_folds = torch.std(tensor, dim=0)
    return mean_folds, sd_folds