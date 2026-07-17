# Supplementary tables

Machine-readable copies of the manuscript supplementary tables that the data
download and analysis workflows depend on. See [`../DATA.md`](../DATA.md) for the
download recipes that use these, and the `hdma` agent skill
([`../.claude/skills/hdma/references/data.md`](../.claude/skills/hdma/references/data.md)).

| file | table | what it is |
|---|---|---|
| `table_s14.tsv` | Table S14 | map of every deposited data type → Zenodo record/DOI (the download index) |
| `table_s2.tsv` | Table S2 | per-cluster cell metadata and the cluster-naming map |
| `table_s6.tsv` | Table S6 | motif lexicon / compendium annotations |

---

## `table_s14.tsv` — Table S14 (data deposition index)

The authoritative map used by the download recipes. One row per deposition; most
data types span several "Part" depositions.

| column | description |
|---|---|
| `data_type` | kind of data (filter on this to select a product) |
| `size` | download size of the deposition |
| `depo_name` | human-readable Zenodo deposition name |
| `record_id` | Zenodo record ID (feed to `get_urls.py`) |
| `url` | concept DOI (resolves to the latest version) |

---

## `table_s2.tsv` — Table S2 (cell metadata / cluster map)

One row per cluster. The `compartment` (col 5) and `Cluster_ChromBPNet` (col 19)
columns drive the cell-type download recipes.

| column | description |
|---|---|
| `Cluster` | unique cluster id |
| `organ` | organ name |
| `organ_code` | abbreviated organ name |
| `cluster_id` | cluster id within organ |
| `compartment` | cellular compartment (epithelial/endothelial/immune/stromal) |
| `L1_annot` | level 1 cell type annotation (n=203) |
| `L2_annot` | level 2 cell type annotation (n=134) |
| `L3_annot` | level 3 cell type annotation (n=37) |
| `dend_order` | order of cluster from hierarchically clustering L1 cell types |
| `ncell` | number of cells in cluster |
| `median_numi` | median RNA number of UMIs per cell |
| `median_ngene` | median RNA number of genes per cell |
| `median_nfrags` | median ATAC number of fragments per cell |
| `median_tss` | median ATAC TSS enrichment ratio per cell |
| `median_frip` | median ATAC fraction of reads in peaks per cell |
| `note` | notes used for cell type annotation |
| `organ_color` | color for clusters based on organ |
| `compartment_color` | color for clusters based on compartment |
| `Cluster_ChromBPNet` | cluster id used in ChromBPNet models |

---

## `table_s6.tsv` — Table S6 (motif lexicon / compendium)

One row per compendium motif. `merged_pattern` maps a compendium motif back to
its per-cell-type constituent motifs.

| column | description |
|---|---|
| `idx_uniq` | unique index for the motif within the compendium |
| `motif_name` | unique motif name (prefix = `idx_uniq`, suffix = `annotation`) |
| `motif_name_safe` | computationally safe `motif_name` (no slashes or pipe) |
| `class` | motif contribution to accessibility (Activating or Repressive) |
| `total_instances_across_celltypes` | total motif instances across all 189 modeled cell types |
| `cwm_logo_fwd` | filename of the forward CWM logo PNG (see Table S14) |
| `cwm_logo_rev` | filename of the reverse CWM logo PNG (see Table S14) |
| `query_consensus` | consensus sequence for the trimmed CWM |
| `annotation` | granular motif annotation (`<family>:<TF>#<index>` or `<family>:<subfamily>/<alt>#<index>`, Vierstra v2.1 convention; composite constituents joined by underscores) |
| `annotation_broad` | broad motif annotation used to group motifs in analyses/figures |
| `category` | motif category (base, base_with_flanks, homocomposite, heterocomposite, partial, repeat, unresolved) |
| `best_match_TF` | TF associated with the best-matching known motif in databases (see Methods) |
| `best_match_motif` | best-matching known motif in databases |
| `best_match_TOMTOM_qval` | TOMTOM q-value for similarity to the best match |
| `best_match_TFs_in_family` | other TFs in the same family as `best_match_TF` (HOCOMOCO) |
| `best_match_TFs_in_Vierstra_archetype` | other TFs whose motifs share the Vierstra v2.1 cluster of the best match |
| `total_seqlets_across_celltypes` | total seqlets across all constituent motifs aggregated into this compendium motif |
| `merged_pattern` | name of the per-cell-type motif cluster corresponding to this compendium motif (maps back to constituents) |
| `n_constituent_motifs` | number of motifs aggregated into this compendium motif |
| `n_constituent_celltypes` | number of cell types contributing at least one constituent motif |
| `constituent_organs` | number of organs contributing at least one constituent motif |
| `constituent_celltypes` | `L2_annot` of the cell types contributing constituent motifs |
