# Paper numbers

A copy of the generated inputs of the paper *Described eukaryotic species in published species-level phylogenies: placement by molecular data, dates and date uncertainty*
(Richter, bioRxiv doi:10.64898/2026.05.31.729127). The manuscript prints no number that is not
in one of these files.

| File | What it holds |
|---|---|
| `headline_numbers.json` | Every reported quantity with its numerator and denominator. |
| `numbers_macros.tex`, `crown_ages_macros.tex`, `table_datasets_macros.tex` | The numbers as the text prints them. |
| `table_which_tree.csv`, `table_rungs_sensitivity.csv`, `table_worked_example.csv` (and `.tex`) | Tables 1 to 3. |
| `fig1_coverage_cascade.csv`, `fig3_c0c1_bars.csv`, `fig2_coverage_scatter.csv`, `crown_ages.csv`, `table_worked_example.csv` | The data of Figures 1 to 5. |
| `worked_example_bd*.csv` | The birth-death fits of the worked example, written by `code/rebuild/worked_example_bd.R`. |
| `supp_*.csv` (and `.tex`), `supp_macros.tex` | The tables of Supplement S1 to S10 and the numbers its text prints. |
| `coverage_per_partition.csv`, `c0_name_resolution.csv`, `opentree_comparison.csv`, `table_datasets.csv`, `table_s*.csv`, `fig_crown_age_points.csv` | Further per-group tables written by the pipeline. |

`python3 scripts/reproduce_paper.py` recomputes these numbers, those of the Supplement included, from the trees and registers of this
repository and compares them with the files here. `python3 code/rebuild/run_worked_example.py`
refits the worked example (R with ape) and compares the result. The figure scripts in
`code/figures/` draw Figures 1 to 5 from these files into `paper/figures/`.
