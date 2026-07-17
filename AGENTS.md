# AGENTS.md — HDMA-public

Orientation for AI coding agents working in this repository. (Claude Code users
get the same guidance automatically via the shipped `hdma` skill — see
`.claude/skills/hdma/`.)

## What this repo is

The public code + data-access companion for the **Human Development Multiomic
Atlas** (HDMA), a single-cell multiome (RNA + ATAC) atlas of human fetal tissues
with per-cell-type ChromBPNet models, a motif lexicon, and enhancer/variant
analyses. It is a **reproducibility companion**, not a turnkey pipeline. A fresh
clone contains **code + rendered HTMLs only**; all data lives on Zenodo. Genome
build throughout is **hg38 / GENCODE release 42**; analyses assume a Slurm HPC.

## Cross-cutting concepts

- **Data is on Zenodo, not in the repo.** Pull it from
  <https://zenodo.org/communities/hdma>; outputs/data are gitignored.
- **Cluster naming has two schemes:** `cluster_id` / `RNA_Clusters` (Table S2,
  "c"-prefixed) vs. **`Cluster_ChromBPNet`** (Table S2 col 19, e.g. `Adrenal_c0`).
  The latter names almost all ChromBPNet/motif/bigwig/training-region files.
- **Supplementary tables are the keys:** Table S14 (data type → URL/DOI/record ID),
  Table S2 (cell metadata / cluster map), Table S6 (motif lexicon). All three ship
  in [`tables/`](tables) (`table_s14.tsv`, `table_s2.tsv`, `table_s6.tsv`; column
  dictionaries in [`tables/README.md`](tables/README.md)).
- **The repo's own docs are authoritative:** `DATA.md`/`DATA.html`, `README.md`,
  `code/03-chrombpnet/README.md`, `MOTIFS.html`.

## Where to go

| If the task is about… | Read |
|---|---|
| **data** — download/find/load, formats, browser tracks, SRA reprocessing | [`.claude/skills/hdma/references/data.md`](.claude/skills/hdma/references/data.md), then [`DATA.md`](DATA.md) |
| **code** — run/adapt a step, envs, figure→code, repo layout, path config | [`.claude/skills/hdma/references/code.md`](.claude/skills/hdma/references/code.md), then [`README.md`](README.md) |

## Three footguns to know up front

1. **Data isn't in the repo** — download from Zenodo via `tables/table_s14.tsv` +
   the `get_urls.py` helper (see the data reference and `DATA.md`).
2. **Create `code/ROOT_DIR.txt`** (gitignored) with the repo's absolute path
   before running any script — scripts read it instead of hardcoding paths.
3. **Cluster files are named by `Cluster_ChromBPNet`** (e.g. `Adrenal_c0`), not
   `cluster_id`.
