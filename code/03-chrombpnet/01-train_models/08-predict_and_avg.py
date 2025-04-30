# Modified from https://github.com/kundajelab/chrombpnet/blob/master/chrombpnet/evaluation/make_bigwigs/predict_to_bigwig.py

import argparse
import pyBigWig
import numpy as np
import pandas as pd
import pyfaidx
from tensorflow.keras.utils import get_custom_objects
from tensorflow.keras.models import load_model
import tensorflow as tf
import chrombpnet.evaluation.make_bigwigs.bigwig_helper as bigwig_helper
import chrombpnet.training.utils.losses as losses
import chrombpnet.training.utils.data_utils as data_utils 
import chrombpnet.training.utils.one_hot as one_hot
import h5py
import json

NARROWPEAK_SCHEMA = ["chr", "start", "end", "1", "2", "3", "4", "5", "6", "summit"]

def write_predictions_h5py(output_prefix, profile, logcts, coords):

    # open h5 file for writing predictions
    output_h5_fname = "{}_predictions.h5".format(output_prefix)
    h5_file = h5py.File(output_h5_fname, "w")
    # create groups
    coord_group = h5_file.create_group("coords")
    pred_group = h5_file.create_group("predictions")

    num_examples=len(coords)

    coords_chrom_dset =  [str(coords[i][0]) for i in range(num_examples)]
    coords_center_dset =  [int(coords[i][1]) for i in range(num_examples)]

    dt = h5py.special_dtype(vlen=str)

    # create the "coords" group datasets
    coords_chrom_dset = coord_group.create_dataset(
        "coords_chrom", data=np.array(coords_chrom_dset, dtype=dt),
        dtype=dt, compression="gzip")
    coords_start_dset = coord_group.create_dataset(
        "coords_center", data=coords_center_dset, dtype=int, compression="gzip")

    # create the "predictions" group datasets
    profs_dset = pred_group.create_dataset(
        "profs",
        data=profile,
        dtype=float, compression="gzip")
    logcounts_dset = pred_group.create_dataset(
        "logcounts", data=logcts,
        dtype=float, compression="gzip")

    # close hdf5 file
    h5_file.close()


# need full paths!
def parse_args():
    parser = argparse.ArgumentParser(description="Make model predictions on given regions and output to bigwig file.Please read all parameter argument requirements. PROVIDE ABSOLUTE PATHS!")
    parser.add_argument("-cm", "--chrombpnet-model", type=str, action="append", required=True, help="Path to chrombpnet model h5")
    parser.add_argument("-r", "--regions", type=str, required=True, help="10 column BED file of length = N which matches f['projected_shap']['seq'].shape[0]. The ith region in the BED file corresponds to ith entry in importance matrix. If start=2nd col, summit=10th col, then the input regions are assumed to be for [start+summit-(inputlen/2):start+summit+(inputlen/2)]. Should not be piped since it is read twice!")
    parser.add_argument("-g", "--genome", type=str, required=True, help="Genome fasta")
    parser.add_argument("-c", "--chrom-sizes", type=str, required=True, help="Chromosome sizes 2 column tab-separated file")
    parser.add_argument("--output-prefix", type=str, required=True, help="Output bigwig file prefix")
    parser.add_argument("--output-key", type=str, required=True, help="Output file key to specify the type of predictions either 'nobias' or 'uncorrected'")
    parser.add_argument("--output-prefix-stats", type=str, default=None, required=False, help="Output stats on bigwig")
    parser.add_argument("-b", "--batch-size", type=int, default=64)
    parser.add_argument("-t", "--tqdm", type=int,default=0, help="Use tqdm. If yes then you need to have it installed.")
    parser.add_argument("-d", "--debug-chr", nargs="+", type=str, default=None, help="Run for specific chromosomes only (e.g. chr1 chr2) for debugging")
    parser.add_argument("-bw", "--bigwig", type=str, default=None, help="If provided .h5 with predictions are output along with calculated metrics considering bigwig as groundtruth.")
    parser.add_argument("-ob", "--output-bed", type=bool, default=False, help="Whether to ouput bed file with consolidated predicted logcounts per region.")
    args = parser.parse_args()
    print(args)
    return args


def softmax(x, temp=1):
    norm_x = x - np.mean(x,axis=1, keepdims=True)
    return np.exp(temp*norm_x)/np.sum(np.exp(temp*norm_x), axis=1, keepdims=True)

def load_model_wrapper(model_hdf5):
    # read .h5 model
    custom_objects={"multinomial_nll":losses.multinomial_nll, "tf": tf}    
    get_custom_objects().update(custom_objects)    
    model=load_model(model_hdf5, compile=False)
    print("got the model")
    model.summary()
    return model

def main(args):

    # load all models
    models = []
    inputlens = []
    outputlens = []
    for model_path in args.chrombpnet_model:
        model = load_model_wrapper(model_hdf5=model_path)
        models.append(model)
        inputlens.append(int(model.input_shape[1]))
        outputlens.append(int(model.output_shape[0][1]))
    
    assert all(inputlen == inputlens[0] for inputlen in inputlens)
    assert all(outputlen == outputlens[0] for outputlen in outputlens)
    inputlen = inputlens[0]
    outputlen = outputlens[0]

    # load data
    regions_df = pd.read_csv(args.regions, sep='\t', names=NARROWPEAK_SCHEMA)
    print(regions_df.head())
    with pyfaidx.Fasta(args.genome) as g:
        seqs, regions_used = bigwig_helper.get_seq(regions_df, g, inputlen)

    gs = bigwig_helper.read_chrom_sizes(args.chrom_sizes)
    regions = bigwig_helper.get_regions(args.regions, outputlen, regions_used) # output regions

    # subset regions
    if args.debug_chr is not None:
        regions_df = regions_df[regions_df['chr'].isin(args.debug_chr)]
        regions = [x for x in regions if x[0]==args.debug_chr]
    regions_df[regions_used].to_csv(args.output_prefix + "_chrombpnet_" + args.output_key + "_preds.bed", sep="\t", header=False, index=False)

    # get predictions
    sum_logits = None
    sum_logcounts = None
    for model in models:
        pred_logits, pred_logcts = model.predict([seqs],
                                            batch_size = args.batch_size,
                                            verbose=True)

        pred_logits = np.squeeze(pred_logits)

        # sum logits and log counts
        if sum_logits is None:
            sum_logits = pred_logits.copy()
            sum_logcounts = pred_logcts.copy()
        else:
            sum_logits = sum_logits + pred_logits
            sum_logcounts = sum_logcounts + pred_logcts

    # average and consolidate across folds
    average_logits = sum_logits / len(models)
    average_prob = softmax(average_logits)

    print(average_logits.shape)
    print(average_prob.shape) 

    average_logcounts = sum_logcounts / len(models)
    average_counts = np.expand_dims(np.exp(average_logcounts)[:,0],axis=1)

    print(average_logcounts.shape)
    print(average_counts.shape) 

    average_pred_profile = average_counts * average_prob

    # output consolidated predicted profile to bigwig
    bigwig_helper.write_bigwig(average_pred_profile,
                            regions,
                            gs,
                            args.output_prefix + "_chrombpnet_" + args.output_key + ".bw",
                            outstats_file=args.output_prefix_stats,
                            debug_chr=args.debug_chr,
                            use_tqdm=args.tqdm)

    # # write h5 file with consolidated predicted counts & profile
    # coordinates = [[r[0], r[-1]] for r in regions]
    # write_predictions_h5py(output_prefix=args.output_prefix,
    #                        profile=average_pred_profile,
    #                        logcts=average_logcounts,
    #                        coords=coordinates)

    # write bed file with logcounts
    if args.output_bed:
        regions_df_used = regions_df[regions_used]
        # add new column at the end
        regions_df_used.insert(10, "logcounts",  average_logcounts)
        regions_df_used.to_csv(args.output_prefix + '_chrombpnet_' + args.output_key + '_preds_w_logcounts.bed', sep="\t", header=False, index=False)


if __name__=="__main__":
    args = parse_args()
    main(args)
