# Figure reproduction code

This directory contains the scripts and auditable inputs used to reproduce the manuscript figures. The current design separates lightweight processed inputs from the heavier tree workspace so the figure code can be audited without reopening large raw tree packages.

---

## Structure

- `fig1_coverage_barplot.R`
- `fig2_crown_age.R`
- `fig3_dated_undated.R`
- `fig4_global_ltt.R`
- `fig5_nested_circles.R`
- processed figure inputs live in `../../data/processed_manuscript_inputs/`

Scripts write their outputs to `../../paper/figures/`.

---

## Dependency model

There are two kinds of figure inputs.

### Fully local processed inputs

Figures 1, 3, and 5 read all of their tabular inputs from `../../data/processed_manuscript_inputs/`.

### Processed manifests plus external trees

Figures 2 and 4 use manifest CSV files stored in `../../data/processed_manuscript_inputs/`, and they read the actual tree files from the clean `../../data/` root through the `TREE_BASE` path configured in the scripts.

This means the manifests are audit-friendly and versionable, but the tree files remain a separate runtime dependency from the processed figure-input tables.

---

## Scripts and expected inputs

| Figure | Script | Main input(s) | External tree dependency | Output |
|--------|--------|---------------|--------------------------|--------|
| 1 | `fig1_coverage_barplot.R` | `../../data/processed_manuscript_inputs/fig1_coverage_data.csv` | no | `../../paper/figures/figure1_barplot_final.png` |
| 2 | `fig2_crown_age.R` | `../../data/processed_manuscript_inputs/fig2_tree_manifest.csv` | yes | `../../paper/figures/crown_age_vs_species.png` |
| 3 | `fig3_dated_undated.R` | `../../data/processed_manuscript_inputs/fig3_dated_undated_counts.csv` | no | `../../paper/figures/figure_dated.png` |
| 4 | `fig4_global_ltt.R` | `../../data/processed_manuscript_inputs/fig4_tree_manifest.csv`, `../../data/processed_manuscript_inputs/fig4_geological_periods.csv` | yes | `../../paper/figures/figure_global_ltt.png` |
| 5 | `fig5_nested_circles.R` | `../../data/processed_manuscript_inputs/fig5_circle_data.csv` | no | `../../paper/figures/figure_nested_circles.png` |

---

## R packages

```r
install.packages(c("ape", "ggplot2", "scales", "ggrepel", "treemapify"))
```

---

## Running the scripts

Run the scripts from this directory.

```bash
Rscript fig1_coverage_barplot.R
Rscript fig2_crown_age.R
Rscript fig3_dated_undated.R
Rscript fig4_global_ltt.R
Rscript fig5_nested_circles.R
```

Figures 2 and 4 resolve `TREE_BASE` from the clean `../../data/` layout and expect the raw tree sources to be present there.

---

## Audit notes

- The processed CSVs in `../../data/processed_manuscript_inputs/` are the auditable figure inputs.
- Figures 2 and 4 still depend on tree files referenced by manifest, not copied into this directory itself.
- If a manuscript number changes, the relevant figure input CSV should usually be updated before editing the plotting code.
- This directory documents the code that consumes the clean data layers; it is not itself a data layer.
