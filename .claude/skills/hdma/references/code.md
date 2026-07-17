# HDMA code route

How to **run, adapt, and navigate** the analysis code. The authoritative
reference is the repo's own [`README.md`](../../../../README.md) (figure→code table,
installation, per-module demo) and, for the deep-learning stack,
[`code/03-chrombpnet/README.md`](../../../../code/03-chrombpnet/README.md). This file
routes you to them and flags the footguns.

**This repo is a reproducibility companion, not a turnkey pipeline.** It enhances
the paper's Methods; there is no single entry point. Numbered modules run in
sequence, and outputs/data are gitignored (a fresh clone is code + rendered HTMLs
only — inputs come from Zenodo, see the data route).

---

## Module map (`code/`, run in order)

- **`utils/`** — shared R helpers + palettes sourced throughout (see below).
- **`01-preprocessing/`** — `01-snakemake/` (per-batch configs for
  bcl → fragments/count matrices), `02-archr_seurat_scripts/` (per-organ
  ArchR/Seurat creation templates) + gitignored `02-preprocessing_per_organ/`
  instantiations, `03-global/` (global peak set, marker genes, peak2gene).
- **`02-global_analysis/`** — global QC/metadata, cell-type dendrogram, TF
  expression (`.Rmd`).
- **`03-chrombpnet/`** — the largest module; has its **own README**.
  `00-prepare_inputs` (prep training data) → `01-train_models` (train/interpret, QC,
  contribs, TF-MoDISco, bigwigs) → `02-compendium` (motif lexicon) →
  `03-syntax` (cooperativity/synergy).
- **`04-enhancers/`** — caCRE export, fragments→tagalign, ABC model,
  co-accessibility, peak-to-gene, ENCODE/VISTA overlaps.
- **`05-misc/`** — `01` create BPCells object, `02` BPCells track plotting,
  `04` ChromBPNet use-cases tutorial.
- **`06-variants/`** — eQTL analysis, g-chromVAR causal-variant analysis,
  ChromBPNet variant scoring + locus plots.

**Naming convention:** numeric prefix + letter suffix encodes execution order
within a module (`04a-`, `04b-`, …). A common pattern is `NN-<task>.sh` (the
command) paired with `NN-jobscript.sh` (the Slurm submission wrapper).

---

## Path config — the critical footgun for cloners

Scripts do **not** hardcode the repo location. They read small **gitignored**
files that a cloner must create themselves before anything runs:

- **`code/ROOT_DIR.txt`** — absolute path to the repo (e.g. Python
  `open("../../ROOT_DIR.txt")`; also used by the `.smk`). **Create this first.**
- `code/SJ_RENV_DIR.txt`, `code/AK_PROJ_DIR.txt`, `code/03-chrombpnet/config.sh`
  — additional path/config files, also gitignored.

R notebooks additionally use `here::here()`. If a script errors immediately on a
missing path file, this is why.

---

## Environments

- **R 4.1.2.** No `renv.lock`/`DESCRIPTION` is committed — exact package versions
  are recorded in the **"Session info"** section of each rendered HTML (the
  substitute for a lockfile).
- **Python: conda**, one env per analysis. The env **name is embedded in each
  Slurm submission script**, not centrally declared. Package lists live in
  [`code/envs/`](../../../../code/envs) as `conda list` output (chrombpnet, finemo,
  gimmemotifs, modiscolite, tangermeme) — these are **package lists, not
  `conda env create`-able YAMLs**; use them to match versions, not to recreate
  directly.
- Everything assumes a **Slurm HPC on Linux** (ran on CentOS 7.9).

---

## How to run

- **R:** open the RStudio project (`HDMA-public.Rproj`) and knit the `.Rmd`, or
  `source()` the `.R` scripts.
- **Python / ChromBPNet:** run the numbered `.sh` scripts, usually via their
  paired `NN-jobscript.sh` Slurm wrappers, inside the conda env named in the
  script.
- Runtime is large (README notes ~several days on HPC even for one organ); some
  analyses are integrative and cannot run on a single organ.

---

## Reusable helpers

- **R** (`code/utils/*.R`): `hdma_palettes.R`, `plotting_config.R` (central
  theming/palettes); `archr_helpers.R`, `seurat_helpers.R`, `matrix_helpers.R`,
  `misc_helpers.R`, `sj_scRNAseq_helpers.R`, `GO_wrappers.R`,
  `chrombpnet_utils.R`; track plotting `track_helpers_SJ.R`, `track_helpers_BL.R`,
  `trackplots.R`.
- **Python** (`code/03-chrombpnet/`): `chrombpnet_utils/` (io, plot, modisco,
  evaluation, utils) and `tangermeme_utils/` (wrappers, syntax_helpers, eval,
  constants, utils — for *in silico* syntax/synergy experiments).

---

## Finding the code behind a figure

`README.md` has a **figure → code** table mapping each manuscript panel to its
rendered-HTML analysis link and repo script/notebook path. Start there:
[`README.md` → Code to produce the figures](../../../../README.md#code-to-produce-the-figures).

---

## Vignettes / entry points for newcomers

- `code/05-misc/04-ChromBPNet_use_cases.ipynb` — load trained models, make
  predictions, score variants (rendered:
  <https://greenleaflab.github.io/HDMA/code/05-misc/04-ChromBPNet_use_cases.html>).
- `code/05-misc/02-bp_cells_plotting_examples.Rmd` — genome-browser-style track
  plotting from the BPCells object.
- Downloading data: see the data route and [`DATA.md`](../../../../DATA.md).
