#!/usr/bin/bash
#SBATCH --job-name="04-peaks"
#SBATCH --time=24:00:00
#SBATCH --output=../../logs/03-chrombpnet/00/%x-%j.out
#SBATCH --partition=akundaje
#SBATCH --cpus-per-task=3
#SBATCH --mem=60G

# Purpose: call peaks on pseudoreplicates and total fragments. The peaks are ranked
# by Macs2 pValue and only the top 300,000 are kept. Then we select peaks 
# called in the total set which overlap with at least one peak in both pseudoreps.

# Adapted from Kundaje Lab / Salil Deshpande.


# SETUP ------------------------------------------------------------------------

# load conda environment
eval "$(conda shell.bash hook)"
conda activate chrombpnet

# source configuration variables
source ../config.sh

# maximum number of parallel processes to run (Default: 4)
input_parallel=3

# path to base of working directory for processed data
input_basedir=$base_dir

# path to .tsv file containing chromosome orders for use with bedtools intersect -g
export input_chromsizes=$chromsizes

# path to .bed.gz file containing blacklisted regions within genome
export input_blacklist=$blacklist

# path to cluster fragments in Kundaje oak space
# because sorting was done by 03c-sort_frag.sh on Kundaje lab cluster
export cluster_frags_dir=$cluster_frags_dir

export peaks_dir=$peaks_dir

# SCRIPT -----------------------------------------------------------------------

callpeak () {
	dataset=${1}
	in_dir="${cluster_frags_dir%/}"
	p1_dir="${in_dir}/pseudorep1"
	p2_dir="${in_dir}/pseudorep2"
	pT_dir="${in_dir}/pseudorepT"
	out_dir=${peaks_dir}
	
	# DEBUG:
	echo $dataset
	echo $in_dir
	echo $p1_dir
	echo $p2_dir
	echo $pT_dir
	echo $out_dir

  # call peaks for each pseudorep
	echo "@@ ${dataset} calling pseudorep1 peaks..." &
	macs2 callpeak -t ${p1_dir}/${dataset}__sorted.tsv -f BED -n ${dataset}__pseudoreplicate1 -g hs --outdir ${out_dir} -p 0.01 --shift -75 --extsize 150 --nomodel -B --SPMR --keep-dup all --call-summits &> ${out_dir}/log__${dataset}__pseudorep1.txt &
	echo "@@ ${dataset} calling pseudorep2 peaks..." &
	macs2 callpeak -t ${p2_dir}/${dataset}__sorted.tsv -f BED -n ${dataset}__pseudoreplicate2 -g hs --outdir ${out_dir} -p 0.01 --shift -75 --extsize 150 --nomodel -B --SPMR --keep-dup all --call-summits &> ${out_dir}/log__${dataset}__pseudorep2.txt &
	echo "@@ ${dataset} calling pseudorepT peaks..." &
	macs2 callpeak -t ${pT_dir}/${dataset}__sorted.tsv -f BED -n ${dataset}__pseudoreplicateT -g hs --outdir ${out_dir} -p 0.01 --shift -75 --extsize 150 --nomodel -B --SPMR --keep-dup all --call-summits &> ${out_dir}/log__${dataset}__pseudorepT.txt &
	wait
	echo "@@ ${dataset} finished waiting"

	p1_in=${out_dir%/}/${dataset}__pseudoreplicate1_peaks.narrowPeak
	p2_in=${out_dir%/}/${dataset}__pseudoreplicate2_peaks.narrowPeak
	pT_in=${out_dir%/}/${dataset}__pseudoreplicateT_peaks.narrowPeak
	p1_out=${out_dir%/}/${dataset}__pseudoreplicate1_peaks_top.narrowPeak
	p2_out=${out_dir%/}/${dataset}__pseudoreplicate2_peaks_top.narrowPeak
	pT_out=${out_dir%/}/${dataset}__pseudoreplicateT_peaks_top.narrowPeak
	
	# max number of peaks
	npeaks=300000
	
	# DEBUG:
	echo $p1_in
	echo $p2_in
	echo $pT_in
	echo $p1_out
	echo $p2_out
	echo $pT_out

  # get the top 300,000 peaks
	echo "@@ ${dataset} getting top peaks..."
    sort -k 8gr,8gr ${p1_in} | head -n ${npeaks} | sort -k 1,1 -k2,2n > ${p1_out}
    sort -k 8gr,8gr ${p2_in} | head -n ${npeaks} | sort -k 1,1 -k2,2n > ${p2_out}
    sort -k 8gr,8gr ${pT_in} | head -n ${npeaks} | sort -k 1,1 -k2,2n > ${pT_out}

  # only keep peaks from pT which overlap with at least one peak in p1 *and* p2
	echo "@@ ${dataset} intersecting peaks"
    min_overlap=0.5
	chr_order="${input_chromsizes}"
	echo $chr_order
    overlap_output=${out_dir}/${dataset}__peaks_overlap.narrowPeak
    bedtools intersect -u -a ${pT_out} -b ${p1_out} -g ${chr_order} -f ${min_overlap} -F ${min_overlap} -e -sorted | bedtools intersect -u -a stdin -b ${p2_out} -g ${chr_order} -f ${min_overlap} -F ${min_overlap} -e -sorted > ${overlap_output}

  # filter out peaks which overlap with blacklist. 
  # N.B.: ideally, one should extend peaks to the 2,114 bp ChromBPNet input field
  # and then filter out peaks where the input field overlaps with any blacklist region.
	echo "@@ ${dataset} filtering blacklist peaks"
	blacklist="${input_blacklist}"
	echo $blacklist
    filtered_output=${out_dir}/${dataset}__peaks_overlap_filtered.narrowPeak
    bedtools intersect -v -a ${overlap_output} -b ${blacklist} > ${filtered_output}

}
export -f callpeak

all_pseudorep_dir="${cluster_frags_dir%/}/pseudorepT"

# for every [cluster]__sorted.tsv file, extract [cluster]
all_sorted=$(ls $all_pseudorep_dir/*__sorted.tsv | xargs -n 1 -I {} basename {} __sorted.tsv)	

datasets=$( for dataset in ${all_sorted[@]}; do
	    
    # check if the done file exists
    # if not, keep the dataset & get list of all unique ones
    done_file="${peaks_dir}/${dataset}__peaks_overlap_filtered.narrowPeak"
    if [[ ! -f "$done_file" ]]; then
        echo "$dataset"
    fi
	
	done | uniq )


echo "@ Calling peaks on cell types: ${datasets}"

parallel --linebuffer -j ${input_parallel} callpeak {} ::: ${datasets}

# DEBUG:
# callpeak Heart_c0

echo "@ done!"
