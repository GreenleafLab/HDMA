---
name: hdma
description: Use for any task involving the HDMA (Human Development Multiomic Atlas) repository — using its data (fragments, count matrices, Seurat/ArchR/BPCells objects, ChromBPNet models, contribution scores, bigwigs, motifs, motif instances, caCREs, ABC loops) or its code (preprocessing, global analysis, ChromBPNet training/interpretation, enhancer/ABC analysis, variant scoring, figure reproduction). Triggers on downloading data from Zenodo, finding files for a given cell type or organ, loading an object, reprocessing raw SRA data, loading tracks in a genome browser, running or adapting a pipeline step, finding which script made a figure, or setting up R/conda envs. This is a router skill — read the relevant file under references/ before acting, because the repo has non-obvious conventions (data lives on Zenodo not in-repo; cluster files are named by Cluster_ChromBPNet; scripts read a gitignored ROOT_DIR.txt).
---

# HDMA

The **Human Development Multiomic Atlas** (HDMA) is a single-cell multiome
(RNA + ATAC) atlas of human fetal tissues with per-cell-type ChromBPNet models,
a motif lexicon, and enhancer/variant analyses. This repository is a
**reproducibility companion** to the paper — it enhances the Methods; it is not a
turnkey pipeline. Genome build throughout is **hg38 / GENCODE release 42**.

This skill is a **router**. Read the relevant reference file **before acting** —
the repo's own docs (`DATA.md`, `README.md`, `code/03-chrombpnet/README.md`,
`MOTIFS.html`) are authoritative; this skill points you to the right part of them
and flags footguns that cause silent mistakes.

## Cross-cutting concepts (read first)

- **Code-only repo; data on Zenodo.** A fresh clone contains code + rendered
  HTMLs only. All data and analysis products live in the Zenodo community
  <https://zenodo.org/communities/hdma>; outputs/data are gitignored.
- **Cluster naming has two schemes.** `cluster_id` / `RNA_Clusters` (Table S2,
  "c"-prefixed) **vs.** `Cluster_ChromBPNet` (Table S2 column 19, e.g.
  `Adrenal_c0`, `Muscle_c6`). The latter names almost every ChromBPNet / motif /
  bigwig / training-region file — use it when locating those.
- **Supplementary tables are the keys.** Table S14 (data type → URL/DOI/record
  ID), Table S2 (cell metadata / cluster map), Table S6 (motif lexicon
  annotations). All three ship in the repo's `tables/` directory
  (`tables/table_s14.tsv`, `table_s2.tsv`, `table_s6.tsv`; column dictionaries in
  `tables/README.md`).
- **Authoritative docs live in the repo:** `DATA.md` / `DATA.html`, `README.md`,
  `code/03-chrombpnet/README.md`, `MOTIFS.html`.

## Task → reference

| If the task is… | Read |
|---|---|
| download/find/load any deposited data; understand a data format; reprocess raw SRA data; load tracks in a genome browser | [references/data.md](references/data.md) |
| run or adapt a pipeline step; set up R/conda envs; find the script behind a figure; use utils/helpers; understand repo layout & path config | [references/code.md](references/code.md) |

## Conventions used throughout

- Genome build **hg38 / GENCODE 42**; fragment coordinates are **0-based**.
- Assumes a **Slurm HPC on Linux** (CentOS 7.9); R 4.1.2; one conda env per
  Python analysis (env name embedded in each Slurm submission script).
- Outputs, figures, motifs, and data directories are **gitignored** — inputs come
  from Zenodo.
