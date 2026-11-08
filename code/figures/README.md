# Figure scripts

The scripts that draw the six main-text figures of the paper from `paper/`,
`registers/coverage_basis.csv` and, for two of them, the site's files in `site/`. Each checks
what it draws against `paper/headline_numbers.json` and stops if they disagree. Output goes to
`paper/figures/`.

| Figure | Script | Data |
|---|---|---|
| 1 The monthly procedure that keeps the catalogue current | `fig_living_catalogue.py` | `paper/headline_numbers.json`, `CITATION.cff` |
| 2 The groups around the tree of life, share named | `fig_tree_of_groups.py` | `registers/coverage_basis.csv`, `site/data/tree_of_groups.json`, `site/app.js` |
| 3 The four levels as the site's dot field, and the estimated total | `fig_dot_field.py` | `registers/coverage_basis.csv`, `paper/fig1_coverage_cascade.csv`, `paper/fig1_coverage_bars.csv`, `site/` |
| 4 Named coverage against the molecular fraction | `fig_coverage_scatter.py` | `paper/fig2_coverage_scatter.csv` |
| 5 Crown ages and archived intervals | `fig_crown_chronograms.py` | `paper/crown_ages.csv` |
| 6 Worked example | `fig_worked_example.py` | `paper/table_worked_example.csv` |

```sh
pip install matplotlib
python3 code/figures/fig_dot_field.py   # and the other five
```

`atlas_figure_common.py` holds the shared palette and the cross-checks; it finds `paper/` and
`registers/` from its own location. Colour has one meaning across the set: brass, species at the
level a figure shows; ink, described species and marks; greys, absent or outside the count;
blue, a step decided by the author; red, the misapplied sampling fraction of the worked example.
