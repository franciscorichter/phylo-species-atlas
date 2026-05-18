---
partition: Crustaceans
category: Arthropods
audit_version: 1
status: in_progress
last_audited: 2026-05-18
auditor: francisco

estimate_source:
  key: Appeltans 2012
  doi: 10.1016/j.cub.2012.09.036
  doi_verified: unchecked
  year_published: 2012
  cited_value: 150000

canonical_tree:
  filename: crustaceans.nwk
  provenance_group: crustaceans
  study: Wolfe JM et al. (2019)
  doi: 10.1098/rspb.2019.0079
  doi_verified: true
  crossref_match: true

methods:
  inference: bayesian
  software: [Bayesian + ML pipeline]
  data_type: UCE
  n_loci: 410
  dating_method: "410-loci anchored hybrid enrichment; 94 species across 58 decapod families."
  support_metric: posterior_probability
  divergence_ci_available: true
  note: "Decapoda-only; not representative of all crustacean diversity. 95 tips for 77k described crustaceans = 0.1% coverage."

shipped_subclades: []
candidate_subclades:
  - taxon: "Whole-tree replacement — broader crustacean phylogeny"
    rationale: "Wolfe 2019 ships under 'crustaceans' but covers only Decapoda (94 of ~9,000 decapod species; ~0.1% of all crustaceans). Atlas needs a Pancrustacea-level tree for adequate partition representation."
    candidate_papers: []
    priority: high
    status: proposed

resolutions:
  - issue: "Shipped tree (Wolfe 2019, 95 tips) is decapod-only — drastically under-represents the Crustacea partition (~77k described species)."
    resolved_by: "Flagged as high-priority replacement candidate. Status kept as in_progress (not verified) until a broader crustacean tree is identified or the partition is restructured."

website:
  surface_methods: true
  surface_uncertainty: true
  surface_candidates: true
  badge: in_progress
---

# Crustaceans — partition audit

**Status: in_progress.** Wolfe 2019 (Proc Roy Soc B, DOI-verified) ships as
the crustacean canonical but covers only Decapoda (95 tips). This is a
**partition-tree mismatch**: a decapod-only tree shouldn't represent the full
~77k-species Crustacea partition.

Appeltans 2012 estimate (100k-250k bilateral range) appears to be paper-
published (from the World Register of Marine Species effort). High-priority
improvement: identify a broader crustacean tree or split the partition.
