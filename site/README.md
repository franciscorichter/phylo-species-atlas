# Species Phylogenies Catalogue — browser site

Static companion site of the release: `index.html`, `app.js`, `style.css`, the fonts in
`fonts/` (Libertinus Serif, OFL), `data/tree_of_groups.json`, and a generated `data.json`.
No framework; d3 from a CDN, the figure drawn on a canvas.

## What it shows

Every described species of eukaryote as dots, one dot per 500 species (4,172 dots for the
2,085,869 species of `data.json`), arranged in a ring around the tree of life of the 62 groups.
The same dots move between seven steps: all described species; the ring and its tree; the four
coverage levels of the paper lit in turn (named in a reference tree, placed with molecular data,
in a chronogram in absolute time, date uncertainty available); and the described species against
the estimated total, as nested discs at the same density. A click on a group gives its card (the
paper's numbers for the group, its reference tree, the tree file and the R call) and opens the
reference tree itself, one dot per species, read from `standardized/trees/`.

- **Numbers.** Every number is read from `data.json`, which `build.py` writes from the paper's
  files (`registers/coverage_basis.csv`, `paper/headline_numbers.json`,
  `paper/fig1_coverage_bars.csv`). The page allocates the dots of each level to the groups by the
  largest-remainder method, so a level's dots sum to its total ÷ 500, and a group's dots at one
  level never exceed its dots at the level before. Bacteria and Archaea (GTDB species clusters)
  are drawn as open dots, outside every eukaryotic total.
- **Order of the groups.** `data/tree_of_groups.json` is the induced subtree of the Open Tree of
  Life synthesis for one taxon per group (api.opentreeoflife.org/v3). The five catch-all groups
  (Other vertebrates, Other Mollusca, Other arachnids, Minor animal phyla, Other phyla) are not
  clades; each is attached at the common ancestor of its members.
- **Estimated total.** The last step sums the low, central and high estimate of each of the 60
  eukaryotic groups in `data.json`; Bacteria and Archaea are left out.
- **Trees.** `site/data/partitions/` holds the per-group notes (`info.yaml`) and copies of some
  other redistributed trees identical to their files in `standardized/trees/`, and three trees
  only the site holds: the Basidiomycota chronogram (Varga et al. 2019), the Xiphosura tree
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
