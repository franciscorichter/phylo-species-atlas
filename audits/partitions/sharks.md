---
partition: Sharks
category: Vertebrates
audit_version: 1
status: verified
last_audited: 2026-05-18
auditor: francisco

estimate_source:
  key: Weigmann 2016
  doi: 10.1111/jfb.12874
  doi_verified: true
  crossref_match: true
  year_published: 2016
  cited_value: 1450
  newer_authoritative_known: true
  newer_candidate:
    doi: null
    year: 2025
    study: Chondrichthyan Tree of Life / Shark References (live)
    cited_value: null
    recommendation: watch
    notes: "Live databases (sharkreferences.com, Chondrichthyan Tree of Life) update species counts continuously; Weigmann 2016 is the canonical published checklist."

canonical_tree:
  filename: sharks.nwk
  provenance_group: sharks
  study: Stein RW, Mull CG, Kuhn TS, Aschliman NC et al. (2018)
  doi: 10.1038/s41559-017-0448-4
  doi_verified: true
  crossref_match: true
  newer_authoritative_known: false
  newer_candidate: null

methods:
  inference: ml
  software: [RAxML, treePL 1.0, TreeAnnotator 1.7.5, Polytomy Resolver]
  data_type: supermatrix
  n_loci: null
  n_genes: null
  total_bp: null
  substitution_model: GTR+CAT (GTRCAT)
  dated: true
  dating_method: |
    treePL rate-smoothing (penalized likelihood) on the RAxML topology;
    7 fossil calibrations distributed across the Chondrichthyan phylogeny;
    500 root-node ages sampled → 500 dated trees per calibration scenario.
  fossil_calibrations: 7
  support_metric: bootstrap
  support_typical: "1,000 bootstrap replicates per RAxML run; low overall bootstrap when M3CPN excluded."
  divergence_ci_available: true
  divergence_ci_note: "500-tree posterior credible set (one per root-node age); atlas ships the MCC consensus."
  posterior_n: 500
  bootstrap_n: 1000
  methods_url: https://www.nature.com/articles/s41559-017-0448-4

shipped_subclades: []

candidate_subclades: []

resolutions:
  - issue: "Tree contains 1,192 species but only 610 have direct molecular data; the remaining 582 were placed via Polytomy Resolver from taxonomic constraints."
    resolved_by: "Documented in methods block. Users interpreting species-level placement uncertainty for the 582 taxon-addition species should consult Stein 2018 §Methods."

website:
  surface_methods: true
  surface_uncertainty: true
  surface_candidates: true
  badge: verified

api_calls:
  - service: crossref
    target: 10.1038/s41559-017-0448-4
    ts: 2026-05-18T00:00Z
    outcome: hit
    cache: pending
  - service: crossref
    target: 10.1111/jfb.12874
    ts: 2026-05-18T00:00Z
    outcome: hit
    cache: pending
---

# Sharks — partition audit

**Status: verified.** Stein 2018 (tree) and Weigmann 2016 (estimate) both
DOI-verified. The interval `1300–1450–1600` is the only one among the
first four audited vertebrate partitions that's plausibly paper-published
(bilateral, not equal to described), though the verbatim quote from
Weigmann 2016 needs to be confirmed when the paywalled PDF is accessed.

## Methods

Stein 2018 uses RAxML (ML) on a 610-species molecular alignment, then
treePL for rate-smoothing with 7 fossil calibrations, then Polytomy
Resolver to add 582 taxa from taxonomic constraints. Output: 500-tree
credible set (one per root-node-age sample), MCC consensus ships in the
atlas. Methodologically distinct from the Bayesian-MCMC vertebrate
partitions (Mammals, Turtles).
