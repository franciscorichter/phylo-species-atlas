---
partition: Birds
category: Vertebrates
audit_version: 1
status: verified
last_audited: 2026-05-17
auditor: francisco

# --- ESTIMATE-SOURCE VERIFICATION ---
estimate_source:
  key: IUCN Red List 2025
  doi: null
  doi_verified: na
  crossref_match: na
  openalex_id: null
  year_published: 2025
  cited_value: 11195
  newer_authoritative_known: false
  newer_candidate: null
  notes: |
    Estimate source is the IUCN Red List of Threatened Species (live database).
    No DOI applies; the URL is https://www.iucnredlist.org/. The 11,195-species
    count is a live snapshot. Birds taxonomy is also tracked by Clements
    (eBird/Cornell) and IOC (International Ornithological Congress). McTavish
    2025 explicitly aligns to Clements; the IUCN figure is reasonable for
    "described species" but differs slightly from Clements 2021 (10,824) and
    2023 (11,017). The 10,824 tip count corresponds exactly to Clements 2021.

# --- CANONICAL TREE VERIFICATION ---
canonical_tree:
  filename: birds.nwk
  provenance_group: birds_mctavish
  study: McTavish EJ, Gerbracht JA, Holder MT, Iliff MJ, Lepage D et al. (2025)
  doi: 10.1073/pnas.2409658122
  doi_verified: true
  crossref_match: true
  newer_authoritative_known: false
  newer_candidate: null
  alternative:
    key: Jetz 2012
    doi: 10.1038/nature11631
    note: |
      Atlas's data_provenance.csv carries an alternative canonical from Jetz et al. 2012
      Nature (birds_jetz, 9,993 tips, BEAST supermatrix with Hackett backbone, 1000-tree
      posterior at VertLife). The atlas's standardized release ships the McTavish 2025
      tree as canonical because it (a) is more recent and updatable, (b) covers more
      species (10,824 vs 9,993), and (c) explicitly notes that 14% of branches differ
      from Jetz 2012 — Jetz's classic BEAST tree remains useful for diversification
      analyses that need a posterior credible set.

# --- METHODS & UNCERTAINTY (McTavish 2025) ---
methods:
  inference: synthesis
  software: [Open Tree synthesis algorithm, Chronosynth (Python), python-opentree, addTaxa (R)]
  data_type: supertree
  n_loci: null
  n_genes: null
  total_bp: null
  substitution_model: null
  dated: true
  dating_method: |
    Chronosynth — summarizes node ages from dated input trees. Expanded from
    Datelife concepts; uses input node-age sets as calibrations and smooths
    over the topology. The dated synthesis tree is built by sampling 100
    random taxon-addition trees per Clements taxonomy year, dating each via
    Chronosynth, then summarizing to a single MCC tree.
  fossil_calibrations: null
  support_metric: source_count
  support_typical: |
    Per-branch support is the number of input studies supporting that branch:
    grey (0 conflicting studies) → dull red (1-4) → bright red (5+ conflicting).
    14% of branches conflict with the majority of input studies. 3,781 branches
    (34% of the tree) have at least one conflicting input study.
  divergence_ci_available: true
  divergence_ci_note: |
    CIs on node dates from the 100-tree taxon-addition sample. The synthesis
    explicitly summarizes uncertainty across input studies. The atlas ships
    only the MCC consensus.
  posterior_n: 100
  bootstrap_n: null
  methods_url: https://www.pnas.org/doi/10.1073/pnas.2409658122
  inputs:
    n_published_studies: 262
    publication_year_range: "1990-2024"
    species_with_phylogenetic_information: 9239
    species_added_by_addTaxa: 1585
    note: |
      9,239 / 10,824 species (85%) have direct phylogenetic information from
      input studies. The remaining 15% are placed by addTaxa using
      Clements/Birds-of-the-World taxonomic constraints. This is unusual
      relative to other partitions in the atlas (typically all tips come
      from sequence data); the audit surfaces this fact as a methodological
      consideration when interpreting species-level placement uncertainty.

# --- SUB-CLADE COVERAGE ---
shipped_subclades:
  - filename: parrots.nwk
    study: Burgio KR & Davis E. (2019)
    taxon: Psittaciformes
    ntips: 413
    methods_inherits_parent: false
    note: |
      Burgio 2019 MRP supertree from 188 source trees. 398 extant species
      + 15 extinct species = 413 tips. Atlas ships all 413; for species-level
      coverage of EXTANT diversity the effective tip count is 398. Future
      atlas-derived pruning could surface this distinction explicitly
      (similar to the turtles canonical-tree swap).

candidate_subclades:
  - taxon: Passeriformes (passerine songbirds)
    rationale: |
      ~6,500 passerine species, ~60% of all birds. The largest clade by far;
      the McTavish 2025 synthesis captures all of them via 262 input studies,
      but a dedicated species-level passerine supermatrix tree would carry
      more uniform support estimates. Oliveros 2019 PNAS is the canonical
      multi-locus passerine tree.
    candidate_papers:
      - doi: 10.1073/pnas.1813206116
        year: 2019
        study: Oliveros et al. — Earth history and the passerine superradiation
        ntips: 4060
        url: https://www.pnas.org/doi/10.1073/pnas.1813206116
    download_target: supplementary
    recipe_file: null
    priority: medium
    status: proposed
    notes: |
      4,060 species (~62% of Passeriformes); BEAST-dated multi-locus tree.
      Complementary to the synthesis: dedicated inference gives uniform
      methodology across the clade, while synthesis aggregates many studies.

  - taxon: Trochilidae (hummingbirds)
    rationale: |
      ~360 hummingbird species. McGuire 2014 is the dedicated multi-locus
      phylogeny. Species-richness is comparable to turtles, with concentrated
      diversification in the Andes; useful for atlas comparisons.
    candidate_papers:
      - doi: 10.1016/j.cub.2014.03.016
        year: 2014
        study: McGuire, Witt, Remsen, Corl, Rabosky, Altshuler, Dudley — Molecular Phylogenetics and the Diversification of Hummingbirds
        ntips: 284
        url: https://www.cell.com/current-biology/fulltext/S0960-9822(14)00366-3
    download_target: dryad
    recipe_file: null
    priority: low
    status: proposed
    notes: |
      Already well-covered by McTavish 2025 synthesis; a dedicated tree would
      add focused dating uncertainty but not species coverage.

  - taxon: Furnariida + Tyrannida (suboscine passerines)
    rationale: |
      Suboscine passerines (~1,300 species). Harvey et al. 2020 Science is
      the recent multi-locus species-level tree.
    candidate_papers:
      - doi: 10.1126/science.aaz6970
        year: 2020
        study: Harvey, Bravo, Claramunt, Cuervo, Derryberry et al. — The evolution of a tropical biodiversity hotspot
        ntips: 1287
        url: https://www.science.org/doi/10.1126/science.aaz6970
    download_target: supplementary
    recipe_file: null
    priority: low
    status: proposed
    notes: "Specialist clade; useful if a passerine-focused audit is later commissioned."

# --- RESOLUTION LOG ---
resolutions:
  - issue: "Burgio 2019 parrots tree has no DOI in data_provenance.csv."
    resolved_by: "DOI 10.1016/j.dib.2019.103882 (Data in Brief) verified via Crossref. info.yaml/sub-phylos/parrots/papers/sources.json carries it; paper-side CSV stays unchanged."
  - issue: "Alternative canonical exists — Jetz 2012 (birds_jetz, 9,993 tips, BEAST + posterior credible set)."
    resolved_by: "Atlas chose McTavish 2025 as canonical based on coverage and updatability. Jetz 2012 documented as alternative in info.yaml.tree.alternative for users who need a posterior credible set."
  - issue: "Parrots tree includes 15 extinct species (413 tips, 398 extant)."
    resolved_by: "Logged for now — extant-only pruning (similar to the turtles species-level derivation) is deferred. Effective extant coverage = 398 species; user-visible until pruning ships."

# --- WEBSITE INTEGRATION ---
website:
  surface_methods: true
  surface_uncertainty: true
  surface_candidates: true
  badge: verified

# --- AUDIT TRAIL ---
api_calls:
  - service: crossref
    target: 10.1073/pnas.2409658122
    ts: 2026-05-17T13:30Z
    cache: pending
    outcome: hit
  - service: crossref
    target: 10.1038/nature11631
    ts: 2026-05-17T13:30Z
    cache: pending
    outcome: hit
  - service: crossref
    target: 10.1016/j.dib.2019.103882
    ts: 2026-05-17T13:30Z
    cache: pending
    outcome: hit
  - service: crossref
    target: 10.1073/pnas.1813206116
    ts: 2026-05-17T13:30Z
    cache: pending
    outcome: hit
  - service: crossref
    target: 10.1126/science.aaz6970
    ts: 2026-05-17T13:30Z
    cache: pending
    outcome: hit
  - service: crossref
    target: 10.1016/j.cub.2014.03.016
    ts: 2026-05-17T13:30Z
    cache: pending
    outcome: hit
---

# Birds — partition audit

## Summary

Status: **verified**. Both DOIs in scope (McTavish 2025, Jetz 2012 alternative)
verify cleanly at Crossref. The atlas ships the McTavish 2025 synthesis tree
(10,824 tips matched to Clements 2021 taxonomy) as canonical, with the
Burgio 2019 parrot supertree as the one shipped sub-clade.

## Methods coverage

McTavish 2025 is a **synthesis** (supertree from 262 published studies),
not a primary Bayesian/ML inference. The methods block reflects this:
no loci/genes/substitution model field; instead the synthesis aggregates
input phylogenies via the Open Tree algorithm and dates them via
Chronosynth. Support is measured by input-study count (not posterior
probability), which is a meaningfully different uncertainty quantification
than the other vertebrate partitions ship.

15% of tips (1,585 species) are placed via `addTaxa` from Clements
taxonomic constraints rather than from sequence data. The audit surfaces
this in `methods.inputs` so users interpreting species-level placement
uncertainty are aware.

## Findings

- ✅ **Burgio 2019 parrots DOI missing from CSV** — found via Crossref:
  `10.1016/j.dib.2019.103882` (Data in Brief, 2019). The atlas's
  data_provenance.csv left the DOI column blank for this row.
  info.yaml/papers/sources.json carries the verified DOI; paper-side
  CSV stays unchanged.
- ✅ **Two competing canonicals** — Jetz 2012 (BEAST + 1000-tree posterior)
  remains a valid alternative; documented as `info.yaml.tree.alternative`.
  Atlas chose McTavish based on updatability and coverage.
- 🟡 **Parrots tree includes extinct species** — 413 tips = 398 extant +
  15 extinct. Extant-only pruning deferred (similar mechanism as the
  turtles species-level derivation; one for Phase E).

## Sub-clade gaps

Three candidates proposed; all low/medium priority because the McTavish
2025 synthesis already incorporates the relevant input studies:

1. Passeriformes — Oliveros 2019 (4,060 species, BEAST multi-locus). Medium.
2. Trochilidae — McGuire 2014 (284 species, dedicated multi-locus). Low.
3. Suboscines — Harvey 2020 Science (1,287 species). Low.

## Status: verified

All resolvable findings have an action in `info.yaml`. Extant-only parrot
pruning is queued for Phase E but does not gate verification.
