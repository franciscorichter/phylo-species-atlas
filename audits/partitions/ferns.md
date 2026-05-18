---
partition: Ferns
category: Plants
audit_version: 1
status: verified
last_audited: 2026-05-18
auditor: francisco

estimate_source:
  key: PPG I 2016
  doi: 10.1111/jse.12229
  doi_verified: unchecked
  year_published: 2016
  cited_value: 13500
  newer_authoritative_known: true
  newer_candidate:
    doi: null
    year: 2025
    study: Fern Tree of Life (FToL, ferntreeoflife.org)
    cited_value: null
    recommendation: watch

canonical_tree:
  filename: ferns.nwk
  provenance_group: ferns
  study: Nitta JH, Schuettpelz E, Ramirez-Barahona S, Iwasaki W (2022)
  doi: 10.3389/fpls.2022.909768
  doi_verified: true
  crossref_match: true

methods:
  inference: ml
  software: [RAxML, IQ-TREE, treePL]
  data_type: multi-locus
  dating_method: "Automated pipeline combining plastomes + plastid regions; continuously updated; published as Fern Tree of Life (FToL)."
  support_metric: bootstrap
  divergence_ci_available: true
  note: "ML and undated versions also published. R package `ftolr` provides programmatic access."

shipped_subclades: []
candidate_subclades: []

resolutions:
  - issue: "PPG I 2016 (Pteridophyte Phylogeny Group I) is the taxonomic authority but doesn't publish a true-diversity range — the 12000-13500-15000 bounds are heuristic round-numbers."
    resolved_by: "Interval classified as fully-heuristic in info.yaml."

website:
  surface_methods: true
  surface_uncertainty: true
  surface_candidates: true
  badge: verified
---

# Ferns — partition audit

**Status: verified.** Nitta 2022 (Frontiers in Plant Science, CC-BY) — the
Fern Tree of Life. Live updates available via ferntreeoflife.org. Interval
is fully-heuristic (PPG I is a taxonomic classifier, not an estimate source).
