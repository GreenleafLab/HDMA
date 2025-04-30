#!/usr/bin/bash
#SBATCH --output=../../logs/03-chrombpnet/02/04/04b-get_tomtom-%j.out
#SBATCH -p akundaje,wjg,biochem,sfgf
#SBATCH -t 01:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=10G

# fail explicitly for any errors, nonexistent variables, etc
# https://gist.github.com/mohanpedala/1e2ff5661761d3abd0385e8223e16425#set--u
# set -euo pipefail

source ../config.sh

# load conda env
eval "$(conda shell.bash hook)"
conda activate modiscolite

module load system
module load cairo

compiled_h5=$modisco_comp_dir/modisco_compiled.h5
out_dir=$modisco_comp_dir

echo $compiled_h5
echo $out_dir

python -u 04b-get_tomtom_matches.py --modisco-h5 $compiled_h5 \
    --out-dir ${out_dir} \
    --meme-db ${vierstra_dir}/all.dbs.meme \
    --verbose True
                    
echo "done."
