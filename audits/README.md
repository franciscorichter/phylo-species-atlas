# Partition audits

Per-partition audit framework for the atlas's 61 partitions in
`data_estimates.csv`. Each partition has a markdown file in
`partitions/<slug>.md` with structured YAML frontmatter recording:

- **Reference verification** — DOI live? Crossref/OpenAlex match? newer
  authoritative estimate worth tracking?
- **Canonical-tree verification** — same questions for the tree paper.
- **Methods & uncertainty** — inference type, software, loci, calibration,
  support metric, divergence CIs.
- **Interval provenance** (per bound: described / low / total / high) —
  paper-published, derived from described count, projection, or atlas
  heuristic. Forces explicit attribution of every interval number.
- **Sub-clade coverage** — what we ship already + candidate sub-trees we
  could add, with citations and download recipes.
- **Audit trail** — Crossref/OpenAlex calls recorded.

## Current coverage (61/61)

- **24 verified** — shipped canonical tree + verified DOIs + interval
  provenance classified.
- **1 in_progress** — Crustaceans (decapod-only Wolfe 2019 ships under
  the whole partition; Bernot 2023 Pancrustacea identified as replacement
  candidate but Dryad data blocked behind Cloudflare).
- **36 deferred** — no shipped canonical tree. Each carries a "why no
  tree" note with a candidate-paper roadmap for Phase E.

### Interval provenance distribution (audited partitions)

| Class | Count | Examples |
|---|---|---|
| **paper-published-range** | 16 | Insects (Stork 2018 2.6M–7.8M), Amphibians (Mora 2011 10k–20k), Sharks (Weigmann 2016), Fungi (Hawksworth 2017 2.2M–3.8M) |
| **partly-heuristic** | 4 | Mammals (Burgin 2018 high is a 2050 description-rate projection), Fish, Bivalves, Bryozoa |
| **fully-heuristic** | 5 | Birds, Squamates, Turtles, Bryophytes, Ferns — checklist sources don't publish true-diversity ranges; high bound is atlas-derived |

## File layout

```
audits/
  README.md                          ← this file
  schema.yaml                        ← machine-readable field schema
  partitions/
    <slug>.md                        ← one per partition
  fetched/                           ← API response cache (gitignored)
    crossref/<doi-slug>.json
    openalex/<work-id>.json
```

## Status workflow

```
stub → in_progress → verified → shipped
                 ↘  deferred (no tree)
```

| Status | Meaning |
|---|---|
| `stub` | Generated from CSVs; nothing verified yet. |
| `in_progress` | DOIs partly verified, methods extracted, but open issues remain. |
| `verified` | DOIs verified, methods complete, interval-provenance classified. Ready to display. |
| `shipped` | All accepted candidate sub-trees downloaded, standardized, added to `metadata.csv`. |
| `deferred` | No shipped canonical tree exists; only estimate provenance audited. |

## Atlas-derived canonical trees

Five partitions ship an **atlas-derived species-level tree** rather than
the paper's raw multi-individual MCC. The derivation collapses the
multi-individual tree to one representative per species (smallest numeric
tip ID per `Genus_species`):

| Partition | Raw MCC | Derived (species-level) | Redundancy removed |
|---|---:|---:|---:|
| Turtles    | 593 tips    | 287 tips    | 51.6% |
| Bryozoa    | 720 tips    | 517 tips    | 28.2% |
| Primates   | 301 tips    | 272 tips    | 9.6%  |
| Seed plants| 356,305 tips| 342,466 tips| 3.9%  |
| Fungi      | 1,672 tips  | 1,602 tips  | 4.2%  |

The original raw MCC is preserved as a locally-hosted "multi-individual"
file in `site/data/partitions/<slug>/uncertainty/`. Reproducible via
`scripts/derive_species_tree.py`. The derivation is **transparent** in
each partition's `info.yaml` (`tree.derivation:` block) and on the audit
page; the chart marks atlas-derived rows with ★ vs ◇ for paper-shipped.

## Schema

The frontmatter schema is in `schema.yaml`. Top-level keys:

```yaml
partition:            # must match data_estimates.csv `group`
category:             # Vertebrates | Plants | Arthropods | Other animals | Microbes & protists | Not yet represented
audit_version: 1
status:               # stub | in_progress | verified | shipped | deferred
last_audited:         # YYYY-MM-DD
auditor:
estimate_source:      # block — DOI, verification, newer-candidate flag, interval_provenance
canonical_tree:       # block — same for the tree paper; reason for deferred
methods:              # block — inference, software, data type, dating, support, CIs
shipped_subclades:    # list — trees already in standardized/
candidate_subclades:  # list — proposed sub-trees with citations + download recipes
resolutions:          # list — issues found and how they were resolved
website:              # block — surface flags + badge
api_calls:            # list — audit trail
```

The new `interval_provenance` block (under `estimate_source`) attributes
each interval bound:

```yaml
interval_provenance:
  described:        { type: database | paper, source, verbatim_quote? }
  estimated_total:  { type: paper-published | derived | atlas-derived | heuristic, ... }
  estimated_low:    { type: paper-published | derived-from-described | heuristic, ... }
  estimated_high:   { type: paper-published | projection | heuristic, ... }
  overall_classification: paper-published-range | partly-heuristic | fully-heuristic | derived
  auditor_note:
```

## Website integration

`site/build.py` reads every `audits/partitions/*.md` frontmatter and merges
it with the verified-facts subset in `site/data/partitions/<slug>/info.yaml`
into `site/data.json`. The web app surfaces:

- **Main dotplot** (`/`): tips ◇/★ + described ◯ + estimated ● with low–high
  whiskers. Whisker styling reflects interval provenance class (solid grey
  = paper-published, faded amber = heuristic). Star marker for atlas-derived
  canonicals.
- **Audit page** (`/audit.html`): per-partition cards with status badge,
  interval-provenance table (with verbatim source quotes when available),
  DOI corrections, resolutions log, shipped + candidate sub-clades.
  Filterable by status and category.
- **Tree detail panel** (per tree): methods + uncertainty folded inline
  into the metadata fields; species-level coverage where redundancy was
  collapsed.

## How to start a new partition audit

```bash
# Copy a similar partition as template
cp partitions/sharks.md partitions/<slug>.md
# Edit the partition: field, DOI, methods, interval provenance

# Verify DOIs
curl https://api.crossref.org/works/<doi>
# Search for sub-clade candidates
curl https://api.openalex.org/works?search=<clade>+phylogeny

# Stage tree + papers
mkdir -p site/data/partitions/<slug>/papers
cp standardized/trees/<file>.nwk site/data/partitions/<slug>/tree.nwk
cp ~/Dropbox/.../papers/<paper>.pdf site/data/partitions/<slug>/papers/tree.pdf
# Write info.yaml + papers/sources.json
```

After Phase B's `scripts/audit_partition.py` ships, this becomes:

```bash
python scripts/audit_partition.py init <slug>
python scripts/audit_partition.py verify <slug>
python scripts/audit_partition.py scan <slug>
```

## Conventions

- **Slugs** are lowercase, space → underscore: `Seed plants` → `seed_plants`.
- **Partition slug must match** `data_estimates.csv` `group` field
  (normalized).
- **No filenames in markdown body** — reference shipped trees by their
  `filename:` field. Cross-check against `standardized/metadata.csv`.
- **Citations live in `data_sources.csv`**, not in audit bodies. Audits
  reference by `source_key`.
- **Paper-side data is read-only** — `standardized/`, `data_*.csv`. All
  corrections live in `site/data/partitions/<slug>/info.yaml`.
- **API responses cached** in `audits/fetched/<service>/<key>.json`
  (gitignored).
