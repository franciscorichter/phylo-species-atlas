# Figure-generation scripts

Scripts that produce the manuscript figures of the Phylogenetic Atlas paper. Deposited from `Paper/code/figures/` of the working manuscript tree so that the dataset repo carries the canonical, citable code for every figure in the paper.

These scripts were written against the working manuscript tree, where the dataset lives in a sister `Data/` directory and figure outputs go to `Paper/figures/`. When running them from inside this repository, adjust input/output paths accordingly (or symlink `Data/` to this repo root). The scripts are deposited verbatim — they have **not** been modified for the repo layout.

## Active scripts

| Manuscript figure | Script | Inputs | Reads trees? | Output |
|---|---|---|---|---|
| **Figure 1** — Coverage dotplot | `fig1_dotplot.py` | `inputs/fig1_coverage_bars.csv` (small summary CSV: partition, n_described, n_in_atlas) | no | `figure_coverage_dotplot.{png,pdf}` |
| **Figure 2** — Crown age vs species richness | `fig2_crown_age.R` | `inputs/fig2_tree_manifest.csv`, `inputs/fig2_crown_age_hpds.csv`; reads dated trees from `standardized/trees/*.nwk` (resolves tip labels through `standardized/dictionary.csv`) | yes | `crown_age_vs_species.{png,pdf}` |
| **Figure 3** — Atlas growth through time | `fig4_atlas_growth.py` | `inputs/figN_atlas_growth.csv` (per-partition historical record; produced by `code/pipeline/rebuild_growth_csv.py`), `inputs/fig4_tree_manifest.csv` | yes | `figure_atlas_growth.{png,pdf}` |
| **Figure S2** — Data landscape (taxonomy x dating x source) | `figS2_data_landscape.py` | `data_estimates.csv` (repo root: per-partition n_tips, n_species, dated flag, source), `standardized/dictionary.csv` | partial (metadata only) | `figure_data_landscape.{png,pdf}` |
| Coverage-gap helper | `fig_coverage_gaps.R` | `data_estimates.csv`, `data_sources.csv` | no | `figure_coverage_gaps.{png,pdf}` (diagnostic; not in main paper) |
| Crown-age HPD extractor | `extract_crown_age_hpds.R` | trees referenced in `inputs/fig2_hpd_inventory.csv` (BEAST MCC `[&height_95%_HPD=...]` annotations, posterior samples, bootstrap replicates) | yes | `inputs/fig2_crown_age_hpds.csv` (sidecar consumed by `fig2_crown_age.R`) |

`extract_crown_age_hpds.R` is a prerequisite of `fig2_crown_age.R`, not a figure in its own right. Run it whenever the underlying source-tree inventory changes.

## Repo-relative inputs

The scripts expect the following inputs to be reachable from the working directory:

- `standardized/dictionary.csv` — ID -> name resolution (this repo, top-level `standardized/` directory)
- `standardized/trees/*.nwk` — standardized tree files (this repo, top-level `standardized/trees/`)
- `data_estimates.csv` — per-partition summary (this repo root)
- `data_sources.csv` — partition source registry (this repo root)
- `inputs/fig*_*.csv` — small per-figure summary CSVs (live alongside the manuscript, **not** in this repo); produced by the manuscript build pipeline or, for `figN_atlas_growth.csv`, by `code/pipeline/rebuild_growth_csv.py`

## Dependencies

R:
```r
install.packages(c("ape", "ggplot2", "scales", "ggrepel"))
```

Python:
```bash
pip install matplotlib pandas pyyaml
```

## Running

These scripts are deposited as a citable record of how the manuscript figures were produced. To re-run them against this repository you need both the small summary CSVs (`inputs/fig*_*.csv`) shipped with the manuscript and the `standardized/` tree files in this repo. Typical order:

```bash
Rscript extract_crown_age_hpds.R   # only when the HPD inventory changes
python fig1_dotplot.py
Rscript fig2_crown_age.R
python fig4_atlas_growth.py
python figS2_data_landscape.py
```
