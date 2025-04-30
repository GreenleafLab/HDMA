#!/usr/bin/bash
#SBATCH --job-name="04b-signal"
#SBATCH --time=08:00:00
#SBATCH --output=../../logs/03-chrombpnet/00/%x-%j.out
#SBATCH --partition=akundaje,wjg,sfgf,biochem
#SBATCH --cpus-per-task=6
#SBATCH --mem=120G

# Purpose: call peaks on each cluster's fragments (no pseudoreplicates),
# in order to generate the fold-change and p-value signal tracks for visualization.
# https://github.com/macs3-project/MACS/wiki/Build-Signal-Track#user-content-Run_MACS2_bdgcmp_to_generate_foldenrichment_and_logLR_track
# https://github.com/macs3-project/MACS/blob/b9b5499c6370fc2019a7d545b3d1fc1668f4587b/docs/bdgcmp.md
# https://docs.google.com/document/d/1f0Cm4vRyDQDu0bMehHD7P7KOMxTOP-HiNoIvL1VcBt8/edit
#
# Note, I installed ucsc-bedclip and ucsc-bedGraphToBigWig inside the chrombpnet
# conda env to make these commands available as required for this script.


# SETUP ------------------------------------------------------------------------

# load conda environment and modules
eval "$(conda shell.bash hook)"
conda activate chrombpnet
module load biology bedtools samtools

# source configuration variables
source ../config.sh

# maximum number of parallel processes to run (Default: 4)
input_parallel=6

# path to base of working directory for processed data
input_basedir=$base_dir

# path to .tsv file containing chromosome orders for use with bedtools -g
export input_chromsizes=$chromsizes

export signal_dir=$bigwigs_signal_dir

export cluster_frags_dir=${cluster_frags_dir}


# SCRIPT -----------------------------------------------------------------------

callpeak () {
  
  # SET INPUTS
  dataset=${1}
  in_dir="${cluster_frags_dir%/}/fragments"
  out_dir=${signal_dir}
  fragments=${in_dir}/${dataset}__sorted.tsv

  fc_bedgraph="${dataset}.fc.signal.bedgraph"
  fc_bedgraph_srt="${dataset}.fc.signal.srt.bedgraph"    
  fc_bigwig="${dataset}.fc.signal.bw"

  pval_bedgraph="${dataset}.pval.signal.bedgraph"
  pval_bedgraph_srt="${dataset}.pval.signal.srt.bedgraph"    
  pval_bigwig="${dataset}.pval.signal.bw"

  # DEBUG:
  echo "@@ cluster: $dataset"
  echo "@@ frags: $fragments"
  echo "@@ out: $out_dir"
  echo "@@ producing..."
  echo $fc_bedgraph
  echo $fc_bedgraph_srt
  echo $fc_bigwig
  echo $pval_bedgraph
  echo $pval_bedgraph_srt
  echo $pval_bigwig

  echo "@@ ${dataset} calling peaks..." &

  # PEAK CALLING ---------------------------------------------------------------

  # call peaks on the cluster's fragments file
  macs2 callpeak -t ${fragments} -f BED -n ${dataset} \
    -g hs --outdir ${out_dir} \
    -p 0.01 --shift -75 --extsize 150 --nomodel -B --SPMR --keep-dup all --call-summits \
    &> ${out_dir}/log__${dataset}.txt


  # the next steps are adapted from the ENCODE bulk ATAC pipeline at
  # https://docs.google.com/document/d/1f0Cm4vRyDQDu0bMehHD7P7KOMxTOP-HiNoIvL1VcBt8/edit

  # FOLD-ENRICHMENT ------------------------------------------------------------

  # run bdgcmp to get the fold change track, then account for chromosome sizes,
  # making sure features do not extend past the chromosomes
  macs2 bdgcmp -t ${out_dir}/${dataset}_treat_pileup.bdg -c ${out_dir}/${dataset}_control_lambda.bdg \
    --o-prefix ${dataset} --outdir ${out_dir} -m FE
  slopBed -i ${out_dir}/${dataset}_FE.bdg -g $input_chromsizes -b 0 \
    | bedClip stdin $input_chromsizes ${out_dir}/${fc_bedgraph}

  # sort and make bigwig
  # sorting as required by bedGraphToBigWig:
  # https://github.com/ENCODE-DCC/kentUtils/blob/master/src/utils/bedGraphToBigWig/bedGraphToBigWig.c#L43
  sort -k1,1 -k2,2n ${out_dir}/${fc_bedgraph} > ${out_dir}/${fc_bedgraph_srt}
  bedGraphToBigWig ${out_dir}/${fc_bedgraph_srt} ${input_chromsizes} ${out_dir}/${fc_bigwig}

  # clean up
  rm -f ${out_dir}/${dataset}_FE.bdg
  rm -f ${out_dir}/${fc_bedgraph} ${out_dir}/${fc_bedgraph_srt}


  # PVAL -----------------------------------------------------------------------

  # sval counts the number of tags per million in the (compressed) BED file
  sval=$(wc -l <(zcat -f ${fragments}) | awk '{printf "%f", $1/1000000}')
  echo $sval

  macs2 bdgcmp -t ${out_dir}/${dataset}_treat_pileup.bdg -c ${out_dir}/${dataset}_control_lambda.bdg \
    --o-prefix ${dataset} --outdir ${out_dir} -m ppois -S ${sval}
  slopBed -i ${out_dir}/${dataset}_ppois.bdg -g $input_chromsizes -b 0 \
    | bedClip stdin $input_chromsizes ${out_dir}/${pval_bedgraph}

  sort -k1,1 -k2,2n ${out_dir}/${pval_bedgraph} > ${out_dir}/${pval_bedgraph_srt}
  bedGraphToBigWig ${out_dir}/${pval_bedgraph_srt} $input_chromsizes ${out_dir}/${pval_bigwig}

  # clean up
  rm -f ${out_dir}/${dataset}_control_lambda.bdg ${out_dir}/${dataset}_treat_lambda.bdg
  rm -f ${out_dir}/${dataset}_ppois.bdg
  rm -f ${out_dir}/${pval_bedgraph} ${out_dir}/${pval_bedgraph_srt}
  rm -f ${out_dir}/${dataset}_peaks.xls
  echo "@@ done ${dataset}."

}
export -f callpeak

all_frags_dir="${cluster_frags_dir%/}/fragments"

# for every [cluster]__sorted.tsv file, extract [cluster]
all_sorted=$(ls $all_frags_dir/*__sorted.tsv | xargs -n 1 -I {} basename {} __sorted.tsv)	

datasets=$( for dataset in ${all_sorted[@]}; do

    # check if the done file exists
    # if not, keep the dataset & get list of all unique ones
    done_file="${signal_dir}/${dataset}.pval.signal.bw"
    if [[ ! -f "$done_file" ]]; then
        echo "$dataset"
    fi

	done | uniq )


echo "@ Calling peaks on cell types: ${datasets}"

parallel --linebuffer -j ${input_parallel} callpeak {} ::: ${datasets}

echo "@ done!"
