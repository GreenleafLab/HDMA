#!/bin/bash

set -eo pipefail

# get inputs
dataset=$1
finemo_param1=$2
finemo_param2=$3
motif_annotation=$4
annotation_drop=$5

# dataset="Eye_c11"
# finemo_param="counts_v0.16_a0.7"

# load conda environment
eval "$(conda shell.bash hook)"
conda activate finemo

source ../config.sh

# set up params
bias_params="Heart_c0_thresh0.4"

# inputs
finemo_out1="${hits_unified_scratch}/${dataset}/${finemo_param1}/"
finemo_out2="${hits_unified_scratch}/${dataset}/${finemo_param2}/"
peaks_file=${chrombpnet_peaks_dir%/}/${dataset}__peaks_bpnet.narrowPeak

# outputs
out="${hits_reconciled_scratch%/}/${dataset}/${finemo_param1}"

# set intermediate & output filenames
rec_bed="${out}/hits.reconciled.bed"
out_bed="${out}/${dataset}__hits_unified.${finemo_param1}.reconciled.bed"
out_gz="${out_bed}.gz"

rec_tsv="${out}/hits_unique.reconciled.tsv"
anno_tsv="${out}/hits_unique.reconciled.annotated.tsv.gz"

# make output directory if it doesn't exist
[[ -d "${out}" ]] || mkdir -p ${out}

echo "@ writing reconciled hits to ${out}"

python 06b-reconcile_hits.py \
    --hits ${finemo_out1}/hits_unique.tsv \
    --hits2 ${finemo_out2}/hits_unique.tsv \
    --out-dir ${out} \
    --merged-motif-anno $motif_annotation \
    --recall ${finemo_out1}/seqlet_report.tsv \
    --recall2 ${finemo_out2}/seqlet_report.tsv \
    --cor-threshold 0.9 \
    --overlap-distance 3 \
    --annotation-drop ${annotation_drop}

# annotate hits based on genomic localization
export R_LIBS_USER=$RENV_SJ

# TODO: switch anno_tsv back to rec_tsv
out="${hits_reconciled_dir%/}/${dataset}/${finemo_param1}/"
out_tsv="${out}/hits_unique.reconciled.annotated.tsv.gz"
Rscript 06b-annotate_hits.R ${peaks_file} ${out_tsv} ${out_tsv} ${out}
    
# index and rename hits
if [[ -f "${rec_bed}" ]]; then

    echo "@ found hits; indexing"

    # sort, compress and index hits
    module load biology samtools bedtools

    bedtools sort -i $rec_bed > $out_bed
    bgzip -c $out_bed > $out_gz
    tabix -p bed $out_gz

    # clean up
    rm $rec_bed $rec_tsv $out_bed

fi

