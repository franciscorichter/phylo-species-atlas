# Phylo-Species Atlas — browser site

Static companion site of the release: `index.html`, `app.js`, `style.css`, the audit log
(`audit.html`, `audit.js`) and a generated `data.json`. No framework; Plotly and d3 from CDNs.

## What it shows

- **Coverage of each group**: the paper's four levels (named in the reference tree, placed
  with molecular data, in a chronogram in absolute time, date uncertainty available) for each of
  the 62 groups, with the described species of each group as denominator, and the totals over
  the described eukaryotic species. Every value is copied from `registers/coverage_basis.csv`
  and `paper/headline_numbers.json` and formatted as the paper prints it (one decimal, two
  below 1). Crown ages and intervals come from `paper/crown_ages.csv`; the estimated total of
  each group from `paper/fig1_coverage_bars.csv`.
- **Trees**: every redistributed tree, read from `standardized/trees/` (viewer and download).
  `site/data/partitions/` holds no copy of a reference tree. It holds the per-group notes
  (`info.yaml`: methods, caveats, the provenance of each estimate, as text), copies of some
  other redistributed trees identical to their files in `standardized/trees/`, and three trees
  only the site shows: the Basidiomycota chronogram (Varga et al. 2019), the Xiphosura tree
  (Obst et al. 2012) and the Primates tree with one tip per species, which
  `scripts/derive_species_tree.py` of the paper's repository derives from
  `standardized/trees/primates.nwk`.

## Build

```bash
python3 site/build.py              # write site/data.json and standardized/names/ (PyYAML)
python3 site/build.py --check      # exit 1 if site/data.json is not what the build writes
python3 -m http.server 8765        # serve the repository root
open http://127.0.0.1:8765/site/
```

`site/build.py` lists its inputs in its docstring. The paper's repository writes `site/` whole
(`scripts/build_release.py`) and then runs `site/build.py`; `scripts/reproduce_paper.py`
compares every number of `site/data.json` with the paper's macros and Table 1, and every tree
file the site serves with its file in `standardized/trees/`.

## Deployment

`.github/workflows/pages.yml` runs `site/build.py` and deploys `site/`, `standardized/trees/`
and `standardized/names/` to GitHub Pages on a push to `main`. Tree names are fetched per tree
from `standardized/names/<file>.json`; `dictionary.csv` is not loaded by the page.
