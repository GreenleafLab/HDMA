# Human Developmental Multiome Atlas

![](img/hdma_logo_small.png)

This repository accompanies the preprint **_Dissecting regulatory syntax in human development with scalable multiomics and deep learning_** (Liu\*, Jessa\*, Kim\*, Ng\*, et al, bioRxiv, 2025).

## Contents

- [Codebase](https://github.com/GreenleafLab/HDMA?tab=readme-ov-file#codebase)
- [Code to produce the figures](https://github.com/GreenleafLab/HDMA?tab=readme-ov-file#code-to-produce-the-figures)
- [Data availability](https://github.com/GreenleafLab/HDMA?tab=readme-ov-file#data-availability)
- [Vignettes](https://github.com/GreenleafLab/HDMA?tab=readme-ov-file#vignettes)
- [Citation](https://github.com/GreenleafLab/HDMA?tab=readme-ov-file#citation)



## Codebase

This repository is meant to enhance the Materials & Methods section by providing code for the custom
analyses in the manuscript, in order to improve reproducibility for the main results.
However, it is not a fully executable workflow.

* `code` --> pipelines, scripts, and analysis notebooks for data processing and analysis
  * [`utils`](https://github.com/GreenleafLab/HDMA/tree/main/code/utils) --> contains .R files with custom functions and palettes used throughout the analysis
  * [`01-preprocessing`](https://github.com/GreenleafLab/HDMA/tree/main/code/01-preprocessing)
    * `01-snakemake` --> config files for processing raw bcl files into fragment files and count matrices
    * `02-archr_seurat_scripts` --> per organ preprocessing scripts to create final Seurat objects and ArchR projects
    * `03-global` --> creating global objects (e.g. global peak set, marker genes)
  * [`02-global_analysis`](https://github.com/GreenleafLab/HDMA/tree/main/code/02-global_analysis)
    * `01` --> global QC and metadata visualizations per organ and per sample
    * `02`, `03` --> construction of dendrogram on cell type similarity
    * `04` --> calculate TF expression levels
  * [`03-chrombpnet`](https://github.com/GreenleafLab/HDMA/tree/main/code/03-chrombpnet)
    * detailed README [here](https://github.com/GreenleafLab/HDMA/blob/main/code/03-chrombpnet/README.md)
    * `00` --> prepare inputs for training ChromBPNet models
    * `01` --> train and interpret ChromBPNet models
    * `02` --> assembly of motif compendium/lexicon
    * `03` --> downstream analysis of ChromBPNet models and motif syntax/synergy
  * [`04-enhancers`](https://github.com/GreenleafLab/HDMA/tree/main/code/04-enhancers)
    * `01` --> export global accessible candidate cis-regulatory elements (acCREs)
    * `02` --> convert fragment files to tagalign for running Activity-By-Contact model (ABC)
    * `03` --> ABC workflow config files
    * `04` --> acCREs co-accessibility analysis
    * `05` --> acCREs peak-to-gene linkage analysis
    * `06` --> acCREs ABC enhancer-to-promoter linkage analysis
    * `07` --> overlap of HDMA acCREs with ENCODE v4 cCREs
    * `08`, `09` --> overlap of HDMA acCREs with VISTA enhancers
  * [`05-misc`](https://github.com/GreenleafLab/HDMA/tree/main/code/05-misc)
    * `01` --> create global BPCells object
    * `02` --> examples for plotting tracks using BPCells
    * `04` --> examples for ChromBPNet use cases, including how to load models and make predictions
  * [`06-variants`](https://github.com/GreenleafLab/HDMA/tree/main/code/06-variants)
    * `00` to `03` --> analysis related to eQTLs
    * `04` to `05` --> causal variant analysis with gchromvar
    * `06` --> variant scoring using ChromBPNet models
    * `07a` to `07c` --> plot variant scoring results


## Code to produce the figures

Code to reproduce analyses is saved in `code`. This table contains pointers to code for the key analyses associated with each figure.
Where possible, the links in the Analysis column lead to rendered HTMLs.

| Figure | Analysis | Path |
| --- | -------- | ---- | 
| Fig 1b, Fig S2b,c | Global QC and metadata | [`code/02-global_analysis/01-global_QC.Rmd`](https://github.com/GreenleafLab/HDMA/tree/main/code/02-global_analysis/01-global_QC.Rmd) |
| Fig 1c | Dendrogram and dotplot | [`code/02-global_analysis/02-dendrogram.Rmd`](https://github.com/GreenleafLab/HDMA/tree/main/code/02-global_analysis/02-dendrogram.Rmd) |
| Fig 1c | ChromVAR heatmap | [`code/02-global_analysis/03-dendrogram_chromvar.Rmd`](https://github.com/GreenleafLab/HDMA/tree/main/code/02-global_analysis/03-dendrogram_chromvar.Rmd) |
| Fig 2a-e, Fig S2f | ABC linking of acCREs | [`code/04-enhanceres/06-abc.Rmd`](https://github.com/GreenleafLab/HDMA/tree/main/code/04-enhanceres/06-abc.Rmd) |
| Fig 2f-g, Fig S3a, Fig S4k | Analysis of VISTA-overlapping enhancers | [`code/04-enhancers/09-overlap_VISTA.Rmd`](https://github.com/GreenleafLab/HDMA/tree/main/code/04-enhancers/09-overlap_VISTA.Rmd) |
| Fig S2d-e | Overlap of acCREs with ENCODE CREs | [`code/04-enhancers/07-overlap_ENCODE_cCREs.Rmd`](https://github.com/GreenleafLab/HDMA/tree/main/code/04-enhancers/07-overlap_ENCODE_cCREs.Rmd) |
| Fig 3b, Fig 6a, Fig S5 | Plotting tracks at select loci | [`code/03-chrombpnet/03-syntax/02-plot_tracks.Rmd`](https://github.com/GreenleafLab/HDMA/tree/main/code/03-chrombpnet/03-syntax/02-plot_tracks.Rmd) |
| Fig 3c, Fig S4a,b,i,j | ChromBPNet QC | [`code/03-chrombpnet/01-train_models/03-model_QC.Rmd`](https://github.com/GreenleafLab/HDMA/tree/main/code/03-chrombpnet/01-train_models/03-model_QC.Rmd) and [`code/03-chrombpnet/01-train_models/03b-plot_correlation.ipynb`](https://github.com/GreenleafLab/HDMA/tree/main/code/03-chrombpnet/01-train_models/03b-plot_correlation.ipynb) |
| Fig 3d-e, Fig 6b,d, Fig S4d-f, Fig S5b | Motif lexicon/compendium | [`code/03-chrombpnet/03-syntax/01-motif_compendium`](https://github.com/GreenleafLab/HDMA/tree/main/code/03-chrombpnet/03-syntax/01-motif_compendium.Rmd) | 
| Fig S4g-h | Visualize motif instances | [`code/03-chrombpnet/03-syntax/03-visualize_hits.ipynb`](https://github.com/GreenleafLab/HDMA/tree/main/code/03-chrombpnet/03-syntax/03-visualize_hits.ipynb) | 
| Fig 4, Fig 5a, Fig S6 | Analysis of motif cooperativity | [`code/03-chrombpnet/03-syntax/04c-plot_cooperativity_results.Rmd`](https://github.com/GreenleafLab/HDMA/tree/main/code/03-chrombpnet/03-syntax/04c-plot_cooperativity_results.Rmd) |
| Fig 5b | Context-specific motif cooperativity | [`code/03-chrombpnet/03-syntax/05b-context_specific_cooperativity.Rmd`](https://github.com/GreenleafLab/HDMA/tree/main/code/03-chrombpnet/03-syntax/05b-context_specific_cooperativity.Rmd) |
| Fig 6f, Fig S7 | eQTL enrichment analysis | [`code/06-variants/03-enrichment_test_collate_results.R`](https://github.com/GreenleafLab/HDMA/tree/main/code/06-variants/03-enrichment_test_collate_results.R) |
| Fig 7b | g-chromVAR analysis | [`code/06-variants/04-gchromvar.R`](https://github.com/GreenleafLab/HDMA/tree/main/code/06-variants/04-gchromvar.R) |
| Fig 7c-d | Plot tracks for variants of interest | [`code/06-variants/07b_rs12740374_muscle_endo_CAD.R`](https://github.com/GreenleafLab/HDMA/tree/main/code/06-variants/07b_rs12740374_muscle_endo_CAD.R) and [`code/06-variants/07c_rs113892147_lung_macrophage_asthma.R`](https://github.com/GreenleafLab/HDMA/tree/main/code/06-variants/07c_rs113892147_lung_macrophage_asthma.R) |
| Fig S8 | Plot tracks for all fetal-only variants | [`code/06-variants/07a_plot_fetal_only_hits_variant_scoring_results.R`](https://github.com/GreenleafLab/HDMA/tree/main/code/06-variants/07a_plot_fetal_only_hits_variant_scoring_results.R) |



## Data availability

All data and analysis products (including fragment files, counts matrices, cell annotations, global acCRE annotations, ChromBPNet models, motif lexicon, motif instances, and genomic tracks) are deposited at https://zenodo.org/communities/hdma. A list of all data types and the corresponding URL and DOI is provided in Table S14 of the manuscript.

We provide a detailed description of the main data types deposited on Zenodo at [`DATA.md`](https://github.com/GreenleafLab/HDMA/blob/main/DATA.md),
along with a demonstration of how to programmatically download files of interest.

All genomic tracks are also hosted online for interactive visualization with the WashU
Genome Browser here at this link:
https://epigenomegateway.wustl.edu/browser2022/?genome=hg38&hub=https://human-dev-multiome-atlas.s3.amazonaws.com/tracks/HDMA_trackhub.json.


## Vignettes

We provide a few notebooks with examples of how to interact with HDMA data,
analysis outputs, and trained models:

- How to download specific files or data for specific cell types from across the Zenodo records: [`DATA.md`](https://github.com/GreenleafLab/HDMA/blob/main/DATA.md#downloading-data-from-zenodo)
- Plotting genomic tracks using BPCells: [`code/05-misc/02-bp_cells_plotting_examples.Rmd`](https://github.com/GreenleafLab/HDMA/blob/main/code/05-misc/02-bp_cells_plotting_examples.Rmd)
- Use cases for ChromBPNet models and outputs, including visualizing predicted accessibility and contribution scores at a region of interest, loading models, making new predictions, and predicting variant effect: [`code/05-misc/04-ChromBPNet_use_cases.ipynb`](https://github.com/GreenleafLab/HDMA/blob/main/code/05-misc/04-ChromBPNet_use_cases.ipynb)


## Citation

If you use this data or code, please cite:

Liu\*, Jessa\*, Kim\*, Ng\*, et al. bioRxiv 2025.

