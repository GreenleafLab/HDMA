# HDMA data route

How to **find, download, and understand** the HDMA data and analysis products.
The authoritative, exhaustive reference is the repo's own
[`DATA.md`](../../../../DATA.md) (rendered:
[`DATA.html`](../../../../DATA.html)) — this file routes you to the right part of it
and flags the non-obvious footguns. Read the linked `DATA.md` section before
writing download or loading code.

All data lives in the Zenodo community <https://zenodo.org/communities/hdma>.
Nothing but code + rendered HTMLs is in the repo itself.

---

## Getting data from Zenodo — the core workflow

The map of everything deposited is **Table S14**, shipped as
[`tables/table_s14.tsv`](../../../../tables/table_s14.tsv) (see
[`tables/README.md`](../../../../tables/README.md) for the column dictionary of
every supplementary table). Columns:

| col | field | use |
|---|---|---|
| 1 | `data_type` | filter to the kind of data you want |
| 2 | `size` | download-size budgeting |
| 3 | `depo_name` | human-readable deposition name |
| 4 | `record_id` | feed to the helper script |
| 5 | `url` | concept DOI (resolves to latest version) |

Most data types are **split across several "Part" depositions**, so the workflow
is: filter Table S14 rows by `data_type` → collect their `record_id`s → resolve
them to direct file URLs with the `get_urls.py` helper (full source in
[`DATA.md` → Helper script](../../../../DATA.md#helper-script)) → `grep` the file
list for the organ/cell type you want → `wget -i`.

```bash
# record IDs live in column 4; pipe them into the helper
grep models tables/table_s14.tsv | cut -f 4 | tr '\n' ' ' | python get_urls.py
#   -> writes urls_record.txt (one latest-version record URL per record)
#      and urls_file.txt (one direct download URL per file, across all records)

grep Brain urls_file.txt > urls_brain.txt   # keep just Brain files
wget -i urls_brain.txt
for i in *.gz; do tar -xvf $i; done          # models/instances/bigwigs arrive as archives
```

`DATA.md` gives four worked recipes — read them there rather than reinventing:
- **ChromBPNet models for an organ** (`grep models` → `grep Brain`)
- **Seurat object for a tissue** (`grep Seurat` → `grep Brain`)
- **all fragment files** (`grep Fragments` → `grep fragments`)
- **bigwigs for a compartment** (`grep Bigwig`, then filter by cell type via
  Table S2 — see the Table S2 note below)

**Footgun — `zenodo_get`:** the alternative `zenodo_get` package needs a
*versioned* DOI, **not** the concept DOI in Table S14 (which resolves to the
latest version). `get_urls.py` handles the resolution for you; prefer it.

The compartment/cell-type recipe additionally filters with `tables/table_s2.tsv`
(e.g. `awk -F'\t' '$5 == "imm" { print $19 }' tables/table_s2.tsv`, where col 5 =
`compartment`, col 19 = `Cluster_ChromBPNet`). The supplementary tables used by
these recipes — `table_s14.tsv`, `table_s2.tsv` (cell metadata), and
`table_s6.tsv` (motif lexicon) — all ship in the repo's
[`tables/`](../../../../tables) directory; see
[`tables/README.md`](../../../../tables/README.md) for their column dictionaries.

---

## What each deposition contains

Table S14 is the source of truth for exact `record_id`s, URLs, and sizes; this is
a description of *what's inside* each `data_type` and how it's organized. Files
are named by cluster/organ using the conventions in the SKILL — in particular
**`Cluster_ChromBPNet`** (e.g. `Adrenal_c0`, `Muscle_c6`), not `cluster_id`.
Genome build throughout is **hg38 / GENCODE release 42**.

### Fragments + count matrices — one deposition per organ/tissue (~13)
Adrenal, Brain, Eye, Heart, Liver, Lung (×2 Parts), Muscle, Skin, Spleen,
StomachEsophagus, Thymus, Thyroid. Per **sample**:
- ATAC: `<sample>.fragments.tsv.gz` + `.tbi` tabix index. **Coordinates are
  0-based.**
- RNA: MEX trio `<sample>.barcodes.tsv.gz` / `.features.tsv.gz` / `.matrix.mtx.gz`.
  Barcode format `CL<lib>_<r1>+<r2>+<r3>` (e.g. `CL131_C01+A01+A01`); features are
  Ensembl gene IDs + names. Detail:
  [`DATA.md` → Fragment files and count matrices](../../../../DATA.md#fragment-files-and-count-matrices).

### Seurat objects — 4 Parts
One **Seurat V4** RDS per organ, `<organ>_RNA_obj_clustered_final.rds`. RNA only
(no ATAC), QC-passing cells. Three assays: `RNA` (counts), `SCT`, and
**`decontX`** (ambient-corrected — **use `decontX` for expression viz**, as in
the paper). Metadata maps to Table S2. Detail:
[`DATA.md` → Seurat objects](../../../../DATA.md#seurat-objects).

### ArchR projects — 4 Parts
One zip per organ, `<organ>_ATAC_obj_clustered_peaks_final_decontx.zip` (ATAC for
QC-passing cells with RNA appended). Five matrices: `GeneExpressionMatrix`,
`GeneScoreMatrix`, `MotifMatrix`, `PeakMatrix`, `TileMatrix`. Note
`ATAC_Clusters`/`ATAC_Clusters_Peak` are ATAC-only and **not used in the
manuscript**. Detail:
[`DATA.md` → ArchR projects](../../../../DATA.md#archr-projects).

### BPCells object — single deposition
One global BPCells object (metadata, fragments, decontX'd + raw RNA, GENCODE 42
transcripts, Peak2Gene + ABC loops, global caCRE peaks, all motif instance hits).
**Footgun:** BPCells is **on-disk** — the RDS stores *paths* to fragments/matrices,
so you must update paths after downloading/moving. Created by
`code/05-misc/01-bp_cells_create_obj.R`; plotting examples in
`code/05-misc/02-bp_cells_plotting_examples.Rmd`. Detail:
[`DATA.md` → BPCells object](../../../../DATA.md#bpcells-object).

### ChromBPNet models, and shared bias model — 2 Parts
One tar per cluster (`<cluster>.gz`); QC-passing models only → **189 cell types**.
The **shared bias model** (trained on Heart_c0 fold_0) is in Part 1. Each tar →
**15 h5 files** `<cluster>__<fold>__<model>.h5`, folds 0–4, model ∈
`bias_model_scaled` / `chrombpnet_nobias` / `chrombpnet`. **Use
`chrombpnet_nobias` for bias-corrected predictions** — it's what all downstream
analysis uses. Tutorial: `code/05-misc/04-ChromBPNet_use_cases.ipynb`. Detail:
[`DATA.md` → Trained ChromBPNet models](../../../../DATA.md#trained-chrombpnet-models).

### ChromBPNet counts mean contribution scores (h5) — 10 Parts
One `<cluster>__average_shaps.counts.h5` per cluster (base-resolution DeepLIFT
counts-head contributions, averaged over folds). Detail:
[`DATA.md` → ChromBPNet mean contribution scores](../../../../DATA.md#chrombpnet-mean-contribution-scores).

### Bigwigs — 8 Parts
One zipped dir per cluster, `<cluster>__bigwigs.gz`, containing observed
p-value signal (`__obs_pval_signal.bw`), bias-corrected predicted accessibility
(`__mean_pred_corrected.bw`), and base-resolution contribution scores
(`__mean_counts_contribs.bw`). **Load contribution bigwigs as `dynseq` tracks**
in WashU/UCSC. Detail:
[`DATA.md` → Bigwig tracks](../../../../DATA.md#bigwig-tracks-for-observed-and-predicted-accessibility-and-contrib-scores).

### Per-cluster TF-MoDISco motifs (h5) — single deposition
Per-cluster `<cluster>__counts_modisco_output.h5` (CWMs/hCWMs/seqlets in
tfmodisco-lite format: `pos_patterns`/`neg_patterns` → `pattern_N` →
`sequence`/`contrib_scores`/`hypothetical_contribs`/`seqlets`) plus
`<cluster>__counts_modisco_report/` (trimmed logos + `motifs.html`), and a
`merged_modisco_patterns_map.tsv` mapping the **6,362** per-cell-type motifs to
the aggregated compendium. Detail:
[`DATA.md` → Motif lexicon and motifs per cell type](../../../../DATA.md#motif-lexicon-and-motifs-per-cell-type).

### Cell metadata, caCREs, training regions, instances, motif lexicon, ABC loops — single combined deposition
The "regulatory elements" deposition bundles several products:
- **ChromBPNet training regions** — `all_training_regions.tar.gz` → `4-peaks/all/`,
  per cluster: `<cluster>__peaks_bpnet.bed` (1,000 bp windows, summit offset = 500)
  and `<cluster>__nonpeaks_fold_<0-4>.bed` (GC-matched backgrounds per fold). Also
  the input peaks for g-chromVAR.
- **Motif lexicon / compendium** — the 6,362 motifs aggregated + QC'd to **508**
  (summarized in Table S6; interactive at
  [`MOTIFS.html`](../../../../MOTIFS.html) / <https://greenleaflab.github.io/HDMA/MOTIFS.html>).
  Files: `motif_compendium.modisco_object.h5` (TF-MoDISco format),
  `motif_compendium.PPM.memedb.txt` and `...trimmed.PPM.memedb.txt` (MEME), plus
  logo images. Match motif names via `merged_pattern` (column R of Table S6).
- **Motif instances** — per cluster `<cluster>__instances.bed.gz` (chr, start,
  end, motif_name, hit_score, strand, pattern_class) and
  `<cluster>__instances.annotated.tsv.gz` (extended per-instance annotations).
- **caCREs, cell metadata, ABC loops** (per-L1-cluster TSVs).

Detail:
[`DATA.md` → ChromBPNet training regions](../../../../DATA.md#chrombpnet-training-regions),
[Motif lexicon](../../../../DATA.md#motif-lexicon-and-motifs-per-cell-type),
[Motif instances](../../../../DATA.md#motif-instances).

### Raw data metadata + library-prep protocol (Supp. Note 1) — single small deposition
Metadata for the raw genomic data deposited to SRA (see SRA reprocessing below),
plus the library preparation protocol.

**Footgun — record IDs:** some concept DOIs quoted in `DATA.md` prose differ from
the version record IDs in `table_s14.tsv` (latest-version resolution). **Treat
`table_s14.tsv` as authoritative** for record IDs used with `get_urls.py`.

---

## Browser tracks (WashU)

Load all cluster tracks by adding this hub in the WashU Epigenome Browser
(Tracks → Remote tracks):

```
https://epigenomegateway.wustl.edu/browser2022/?genome=hg38&hub=https://human-dev-multiome-atlas.s3.amazonaws.com/tracks/HDMA_trackhub.json
```

Per-cluster tracks: `signal_pval`, `predictions`, `predictions (uncorrected)`,
`peaks`, `deepshap_counts`, `hits_unified` (Fi-NeMo motif instances),
`abc_loops`. Step-by-step + screenshots:
[`DATA.md` → Genomic tracks on the WashU Genome Browser](../../../../DATA.md#genomic-tracks-on-the-washu-genome-browser).

---

## Reprocessing raw data from SRA

Raw data is deposited **anonymized** (read sequences replaced with reference via
BAMboozle; barcodes preserved) under BioProject
[`PRJNA1402391`](https://www.ncbi.nlm.nih.gov/bioproject/1402391). Reprocess
anonymized FASTQs → per-sample fragments + count matrices with the
**`process-anonymize` branch** of
[`GreenleafLab/shareseq-pipeline`](https://github.com/GreenleafLab/shareseq-pipeline/tree/process-anonymize),
which chains `prep_sra.smk` → `ingest_anonymized.smk` → `shareseq.smk`. Output
cell barcodes (`CL{N}_BC1+BC2+BC3`) match the published Seurat/ArchR metadata, so
you can join by barcode string.

Demultiplexing parameter tables live in
`code/01-preprocessing/01-snakemake/` (`sample_barcodes.tsv`,
`sublibrary_indices.tsv`, per-batch `b{N}_<organ>.yaml`).

**Footgun:** the `Organ` column in those tables is a **processing batch**, not a
tissue — most map 1:1 but Lung was split (`Lung_b12`, `Lung_b18`). Filter on the
**exact** `Organ` value; substring matching collides (`grep Lung` hits both;
`grep Spleen` hits `AdrenalThyroidSpleen` and `SpleenThymus`). Full pipeline
overview, required inputs, and expected outputs:
[`DATA.md` → Processing raw anonymized data from SRA](../../../../DATA.md#processing-raw-anonymized-data-from-sra).
