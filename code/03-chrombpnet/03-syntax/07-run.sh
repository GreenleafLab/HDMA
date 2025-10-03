#!/usr/bin/bash

source ../config.sh

# find which cell types to keep
clusters=( $(awk '{print $1}' ${chrombpnet_models_keep2}) )

# divide clusters into groups of 15
for ((i=0; i<${#clusters[@]}; i+=15)); do
    group=( "${clusters[@]:i:15}" )
    sbatch -J cooc_group_$((i/15)) 07-jobscript.sh ${group[@]}
    sleep 1s
done
