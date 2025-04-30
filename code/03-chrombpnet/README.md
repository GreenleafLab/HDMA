
Overview of the steps involved in training ChromBPNet models on the atlas
and analyzing outputs.

We show how to load and use trained ChromBPNet models in the tutorial at [`code/05-misc/04-ChromBPNet_use_cases.ipnyb`](https://github.com/GreenleafLab/HDMA/blob/main/code/05-misc/04-ChromBPNet_use_cases.ipynb) ([html](https://greenleaflab.github.io/HDMA/code/05-misc/04-ChromBPNet_use_cases.html)).

## Tools used

- ChromBPNet: https://github.com/kundajelab/chrombpnet
- TF-MoDISco (tfmodiscolite implementation): https://github.com/jmschrei/tfmodisco-lite
- gimmemotifs: https://gimmemotifs.readthedocs.io/en/master/
- Fi-NeMo: https://github.com/austintwang/finemo_gpu
- tangermeme: https://github.com/jmschrei/tangermeme


## `chrombpnet_utils`

Python utilities for downstream ChromBPNet related analysis, plotting, etc.


## `tangermeme_utils`

Python utilities for downstream ChromBPNet related analysis, particularly
for interfacing with the [tangermeme](https://tangermeme.readthedocs.io/en/latest/index.html) package.
These are mainly used for the *in silico* experimentation.


## 0. Preparing inputs, `00-inputs`

These scripts pertain to generating the necessary inputs for training ChromBPNet models.

- `01`: get fragments per cluster per sample from ArchR projects
- `02`: process fragments for calling peaks (generating pseudoreplicates)
- `03`: group fragments across samples per cell type, and sort fragments
- `04`: call peaks with macs2 on pseudoreplicates
- `04b`: call peaks on the fragments of each cluster (no pseudoreps), to generate p-value and fold-change signal bigwig tracks
- `05`: format peaks for ChromBPNet
- `06`: make bigwigs
- `07`: make chromosome splits for test/train/validation folds
- `08`: get GC-matched negative peaks
- `09`: tabix-index peaks for viewing in browser
- `10`: randomly subset backgrounds to 100 regions per cell type to use for downstream analysis




## 1. Training models, `01-train_models`

These scripts run the ChromBPNet workflow, which involves training a bias model,
using the bias model to train chromatin accessibility models in each cell type,
calculate contribution scores, generate model predictions, and perform motif discovery.

- `01`: train bias models
- `02`: train ChromBPNet models with a given bias model, initial interpretation with 30k subsample peaks
- `03`: perform model QC and make figures with performance metrics
- `04`: get contribution scores (_aka_ DeepLIFT or DeepSHAP scores) as h5 files per fold per cell type
- `05`: average contribs across folds, separately for profile and counts
- `06`: deep TF-MoDISco run with 1M seqlets on averaged contribs
- `07`: generate contribution score bigwigs, for viewing in browser
- `08`: predict accessibility in peak regions across folds, and average the predictions (bias-corrected) to generate bw
- `09`: prepare supplementary table with ChromBPNet performance and QC metrics





## 2. Assembly of motif compendium, `02-compendium`

These scripts pertain to assembling a motif compendium from across the atlas,
calling instances of each motif in each cell type, and then annotating those instances.

- `01`: convert MoDISco CWMs to PFMs to prepare the input for `gimme cluster`:
this is run _per organ_, and motifs from all cell types (models) within an organ are concatenated.
Done separately for positive and negative patterns.
- `02`: for each organ, run `gimme cluster` on PFMs of motifs from all cell types in that organ, to get broad clusterings.
Done separately for positive and negative patterns.
- `02b`: run `gimme cluster` again, using as input all the clusters from the _per organ_ runs of `gimme cluster` in Step 02.
Dine separately for positive and negative patterns. This produces "superclusters" of patterns from across cell types.
- `03`: merge similar MoDISco patterns within each supercluster together, using the `modiscolite.SimilarPatternCollapser` functionality.
This produces one MoDISco h5 object per supercluster, containing merged, non-redundant patterns.
- `04`: combine MoDISco h5 objects from the merging step into one h5 object, match to known patterns, assign unique index, and prepare as MEME format
- `05`: call hits in each cell type, using the unified motifs
- `06`: generate HTML reports to assess the quality of hit calling with unified motif set, filter & label hits using manual annotation
- `06b`: filter out low QC hits, reconcile overlapping hits, and get genomic annotations
- `07`: run NucleoATAC per cell type


## 3. Analysis of TF binding site syntax, `03-syntax`

- `00`: compute the mean and summed importance score for all trimmed CWMs
- `01`: global overview of de novo motifs and predictive instances
- `02`: visualization of genomic tracks at example loci
- `03`: visualization of motif instances (hits) for a few examples
- `04`: in silico analysis of cooperativity
  - `04a`: assemble the set of composites to test
  - `04b`: run the in silico analysis
  - `04c`: summarize the results and generate plots
  - `04d`: make animations showing hard and soft synergy
- `05`: in silico analysis of cell type-specific cooperativity
  - `05a`: run the analysis
  - `05b`: visualize results


_**NOTE**_: variant scoring using ChromBPNet models is done in the `code/06-variants` directory.



