# Web-app data folder

This folder holds the **curated data the web app serves**. It is distinct
from the paper's data:

| Folder | Role | Edit policy |
|---|---|---|
| `standardized/` (repo root) | The paper's published dataset — 264 Newick trees + dictionaries. | **Do not modify** for any web-app reason. Changes here are paper changes. |
| `data_estimates.csv`, `data_provenance.csv`, `data_sources.csv` (repo root) | The paper's source tables. | **Do not modify** for any web-app reason. |
| `audits/` (repo root) | Working area for partition verification (DOI checks, methods extraction, candidate sub-clade scans). | Edit freely; this is the process. |
| `site/data/` (this folder) | What the web app serves. Self-contained per-partition bundle. | Edit freely; this is the product. |

The pipeline is **paper data + audit → curate → `site/data/` → web app**. When
an audit advances to `verified`, its facts (and the corresponding tree + papers)
get materialised here so the web app can consume them without reaching back
into the paper's working tree.

## Layout

```
site/data/
├── README.md                           ← this file
└── partitions/
    └── <partition-slug>/
        ├── tree.nwk                    ← copy of standardized/trees/<file>.nwk
        ├── info.yaml                   ← partition metadata for the web app
        ├── papers/
        │   ├── tree.pdf                ← the paper that provides the tree
        │   ├── estimate.pdf            ← the paper that provides the diversity estimate
        │   └── sources.json            ← DOI, citation, license, open-access flag
        └── sub-phylos/
            └── <subclade-slug>/
                ├── tree.nwk
                ├── info.yaml
                └── papers/
                    ├── tree.pdf
                    └── sources.json
```

- **`tree.nwk`** — verbatim copy of the standardized tree file. Same numeric
  tip IDs that resolve against `standardized/dictionary.csv`.
- **`info.yaml`** — the verified subset of the partition audit (category,
  estimate value + range + source key, canonical-tree provenance, methods
  block). The web app reads this directly.
- **`papers/`** — the two PDFs (one for the tree, one for the estimate),
  plus `sources.json` mapping each file to its DOI/citation/license.
- **`sub-phylos/<slug>/`** — one folder per sub-clade tree, same shape as
  the partition folder (tree + info + papers).

## PDF policy

`papers/*.pdf` is **gitignored by default** — see `site/data/.gitignore`.
PDFs exist locally and can be included in deploys, but are not committed
to GitHub unless explicitly excepted.

To commit a specific PDF (e.g. a CC-BY paper from PLOS / MDPI / eLife), add
a negation rule to `site/data/.gitignore`:

```
!partitions/mammals/papers/tree.pdf
```

`sources.json` records the licensing status (`open_access: true/false`)
so the build script can decide which PDFs to bundle into a public deploy.

## How to add a new partition

1. Run the audit (`audits/partitions/<slug>.md`) and advance to `verified`.
2. Create `site/data/partitions/<slug>/`.
3. Copy the canonical Newick file from `standardized/trees/` to `tree.nwk`.
4. Copy the two papers from the local PDFs collection into `papers/`,
   named `tree.pdf` and `estimate.pdf`.
5. Write `info.yaml` from the audit's verified frontmatter.
6. Write `papers/sources.json` with DOI/citation/license for each paper.
7. For each shipped sub-clade tree, repeat steps 3–6 inside
   `sub-phylos/<subclade-slug>/`.

After Phase B (`scripts/audit_partition.py`) ships, steps 2–7 become a
single command:

```bash
python scripts/audit_partition.py materialize <slug>
```

…which reads the verified audit, copies the tree(s), and writes the
`info.yaml` + `sources.json` files.

## Web app integration

`site/build.py` (extended in Phase D) reads `site/data/partitions/*/info.yaml`
to assemble `data.json` instead of joining the paper CSVs directly. Until
Phase D, the build still uses the paper CSVs; this folder is staging.
