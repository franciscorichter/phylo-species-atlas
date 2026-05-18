---
partition: Bryophytes
category: Plants
audit_version: 1
status: verified
last_audited: 2026-05-18
auditor: francisco

estimate_source:
  key: CoL 2025
  doi: null
  doi_verified: na
  year_published: 2025
  cited_value: 24000
  newer_authoritative_known: false

canonical_tree:
  filename: bryophytes.nwk
  provenance_group: bryophytes
  study: Bechteler J et al. (2023)
  doi: 10.1002/ajb2.16249
  doi_verified: true
  crossref_match: true

methods:
  inference: ml
  software: [ASTRAL, RAxML]
  data_type: UCE
  n_genes: 228
  total_bp: null
  dating_method: "Phylogenomic time tree from 405 exons (228 nuclear genes), 52 of 54 bryophyte orders. 29 fossil calibrations."
  fossil_calibrations: 29
  support_metric: bootstrap
  divergence_ci_available: true
  note: "Also ASTRAL coalescent tree and RAxML tree available."

shipped_subclades: []
candidate_subclades: []

resolutions:
  - issue: "Tree has only 533 tips for ~24,000 described bryophyte species — 2.2% coverage. Severe sampling sparsity."
    resolved_by: "Coverage gap is fundamental to bryophyte phylogenetics (most species have no sequence data). Documented; no fix available without major new sequencing effort."

website:
  surface_methods: true
  surface_uncertainty: true
  surface_candidates: true
  badge: verified
---

# Bryophytes — partition audit

**Status: verified.** Bechteler 2023 (American Journal of Botany) — phylogenomic
time tree using 228 nuclear genes from 405 exons; 29 fossil calibrations; 52/54
bryophyte orders. 533 tips covering 2.2% of ~24k described species — the
highest-quality bryophyte tree available but inherently sparse at species level.

Interval is fully-heuristic (CoL is a checklist database).
