---
partition: Bryozoa
category: Other animals
audit_version: 1
status: verified
last_audited: 2026-05-18
auditor: francisco

estimate_source:
  key: Appeltans 2012
  doi: 10.1016/j.cub.2012.09.036
  doi_verified: true
  year_published: 2012
  cited_value: 8000

canonical_tree:
  filename: bryozoa.nwk
  provenance_group: bryozoa
  study: Orr RJS et al. (2022)
  doi: 10.1126/sciadv.abm7452
  doi_verified: true

methods:
  inference: ml
  data_type: UCE
  note: "17-gene supermatrix, RAxML ML, rogue-taxon pruning; cheilostome-focused."

resolutions:
  - issue: "Orr 2022 MCC ships with 720 tips but only 517 species (27.8% subspecies/specimen redundancy — highest of any tree in the atlas)."
    resolved_by: "Atlas-derived 517-tip species-level pruning is now the canonical tree (one tip per species, smallest numeric ID per Genus_species). The original 720-tip MCC is preserved as a locally-hosted uncertainty file. Reproducible via scripts/derive_species_tree.py."

website:
  badge: verified
---

# Bryozoa — partition audit

**Status: verified.** Orr RJS et al. (2022) tree (DOI 10.1126/sciadv.abm7452). Appeltans 2012 estimate. Interval
classification: partly-heuristic.

Canonical is now the **atlas-derived 517-tip species-level pruning** of
the Orr 2022 MCC (720 → 517, removing 27.8% redundancy from
subspecies/specimen duplicates — the highest redundancy in the atlas).
Coverage 517 / 6000 = 8.6% at species level. Original multi-individual
MCC preserved as uncertainty data.
