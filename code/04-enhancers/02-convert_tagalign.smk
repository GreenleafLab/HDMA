import glob
import os

with open("../ROOT_DIR.txt", "r") as file:
    base_dir = file.readline().strip()
frag_folder = f"{base_dir}/output/03-chrombpnet/00-inputs/cluster_fragments/fragments"
out_folder = f"{base_dir}/output/04-enhancers/02/tagalign/"

frag_files = glob.glob(frag_folder + "/*__sorted.tsv")
all_clusters = [os.path.basename(file).replace("__sorted.tsv", "") for file in frag_files]
all_clusters = sorted(all_clusters)


# test code
#all_clusters = ["Adrenal_c9", "Thyroid_c9"]

rule all:
  input: expand(out_folder + "/{cluster}_tagAlign.gz", cluster=all_clusters)
  output: "success.txt"
  shell:
    "echo \"success\" > {output}"


rule convert_tagalign:
  input:
    frag = frag_folder + "/{cluster}__sorted.tsv"
  output:
    tagalign = out_folder + "/{cluster}_tagAlign.gz",
    tagalign_tbx = out_folder + "/{cluster}_tagAlign.gz.tbi"
  threads: 4
  shell:
    "LC_ALL=C cat {input.frag} | sed '/^#/d' | awk -v OFS='\t' '{{mid=int(($2+$3)/2); print $1,$2,mid,\"N\",1000,\"+\"; print $1,mid+1,$3,\"N\",1000,\"-\"}}' | sort -k 1,1V -k 2,2n -k3,3n --parallel {threads} | bgzip -c > {output.tagalign} ; "
    "tabix -p bed {output.tagalign}"
