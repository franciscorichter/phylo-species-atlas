# Data-standardization pipeline

Pipeline scripts deposited from `Data/scripts/` of the working manuscript tree. These scripts are part of the dataset construction layer (raw → curated → standardized → derived → processed). They are deposited verbatim — paths inside the scripts assume the manuscript layout (a sister `Data/` directory with `standardized/`, `site/`, etc.) and have **not** been rewritten for the dataset repo.

## Deposited scripts

### `derive_species_tree.py`

**Stage:** standardized → site (per-partition representative species tree)

**Purpose:** Given a standardized Newick tree where tips are numeric IDs (resolved through `standardized/dictionary.csv` to specimen-level names of the form `Genus_species_specimen_code`), prune to **one representative tip per species** (species = first two underscore-separated tokens of the standardized name).

**Branch-length policy:** When collapsing siblings of the same species, the representative inherits its own original branch length. Internal unary nodes are collapsed by summing branch lengths into the surviving child.

**Inputs:**
- `--in` standardized Newick tree (e.g. `standardized/trees/<slug>.nwk`)
- `--dictionary` `standardized/dictionary.csv` (ID → name lookup)

**Outputs:**
- `--out` per-partition species-level Newick (e.g. `site/data/partitions/<slug>/tree.nwk`)

**Usage:**
```bash
python derive_species_tree.py \
    --in   standardized/trees/turtles.nwk \
    --dictionary standardized/dictionary.csv \
    --out  site/data/partitions/turtles/tree.nwk
```

### `build_experts.py`

**Stage:** documentation scaffolding (internal expert-agent system; not part of the dataset itself)

**Purpose:** Scaffold per-partition and per-category expert-agent folders under `Data/experts/`. For each partition / category produces seven artifacts: `AGENT.md`, `report.html`, `timeline.csv`, `subclades.yaml`, `discoveries.yaml`, `gaps.md`, `references.bib`, `state.yaml`. Used internally for research-protocol tracking — **not coupled to the parent paper or the live website**, so it is deposited here for completeness only.

**Inputs:** templates under `Data/experts/_shared/`, partition list from `Data/site/data/partitions/`.

**Outputs:** `Data/experts/<slug>/` folder per partition.

**Usage:**
```bash
cd Data
python3 scripts/build_experts.py [--force]
```

### `rebuild_growth_csv.py`

**Stage:** processed (manuscript input)

**Purpose:** Rebuild `figN_atlas_growth.csv` (the historical record consumed by Figure 3) using an **honest historical-counting** rule: a partition appears in year X if any species-level tree meeting atlas inclusion criteria was published by then, not only if today's canonical tree was published by then. Per-partition history is curated from audit MDs, `data_provenance.csv`, and DOI-verified literature.

**Inputs:** in-source `HISTORY` table; reads `data_provenance.csv` and partition audit MDs for cross-checks.

**Outputs:** `Paper/inputs/figN_atlas_growth.csv` (consumed by `code/figures/fig4_atlas_growth.py`).

## Pipeline coverage — what is and is not here

The three deposited scripts cover the **derived** (species-tree extraction) and **processed** (growth CSV) stages, plus an internal documentation scaffolder.

The earlier dataset-construction stages of the pipeline — **tip-label cleaning, integer-ID assignment, Catalogue of Life (CoL) mapping, and the dictionary build** that produces `standardized/dictionary.csv` — are not present as scripts in `Data/scripts/`. They were executed as a sequence of ad-hoc per-partition operations on each source tree before standardization and are not currently captured in a single reproducible pipeline script.

To make the dataset fully reproducible, these earlier stages need to be reconstructed and deposited here as separate scripts. See the audit's items-requiring-manual-action list.
