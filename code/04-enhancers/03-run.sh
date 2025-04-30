# normal
#snakemake -j1

# custom peak set
# complete scripts at https://github.com/GreenleafLab/ABC-Enhancer-Gene-Prediction-CustomRegions/tree/main
snakemake --profile profile -s workflow/Snakefile_custom_peaks.smk

