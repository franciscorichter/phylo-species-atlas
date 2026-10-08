# Figure scripts

The scripts that draw the five main-text figures of the paper from `paper/` and
`registers/coverage_basis.csv`. Each checks what it draws against
`paper/headline_numbers.json` and stops if they disagree. Output goes to `paper/figures/`.

| Figure | Script | Data |
|---|---|---|
| 1 The four coverage levels and the losses between them | `fig_coverage_cascade.py` | `paper/fig1_coverage_cascade.csv` |
| 2 Named and molecular coverage per group | `fig_c0c1_bars.py` | `paper/fig3_c0c1_bars.csv` |
| 3 Named coverage against the molecular fraction | `fig_coverage_scatter.py` | `paper/fig2_coverage_scatter.csv` |
| 4 Crown ages and archived intervals | `fig_crown_chronograms.py` | `paper/crown_ages.csv` |
| 5 Worked example | `fig_worked_example.py` | `paper/table_worked_example.csv` |

```sh
pip install matplotlib
python3 code/figures/fig_coverage_cascade.py   # and the other four
```

`atlas_figure_common.py` holds the shared palette and the cross-checks; it finds `paper/` and
`registers/` from its own location.
