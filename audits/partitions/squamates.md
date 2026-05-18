---
partition: Squamates
category: Vertebrates
audit_version: 1
status: verified
last_audited: 2026-05-18
auditor: francisco

estimate_source:
  key: Uetz 2025
  doi: null
  doi_verified: na
  crossref_match: na
  year_published: 2025
  cited_value: 11744
  newer_authoritative_known: false
  newer_candidate: null
  notes: "Reptile Database (reptile-database.org) is a live taxonomic database. No DOI; Uetz P maintains the canonical reptile species list."

canonical_tree:
  filename: squamates.nwk
  provenance_group: squamates
  study: Tonini JFR, Beard KH, Ferreira RB, Jetz W, Pyron RA (2016)
  doi: 10.1016/j.biocon.2016.03.039
  doi_verified: true
  crossref_match: true
  newer_authoritative_known: false
  newer_candidate: null

methods:
  inference: bayesian
  software: [BEAST (via Pyron 2013), Jetz et al. 2012 imputation method]
  data_type: supermatrix
  n_loci: 12
  data_source_note: |
    Built on Pyron 2013 12-gene supermatrix (4,170 squamates + 1 Sphenodon).
    Tonini 2016 extends to 9,755 species via taxonomic imputation (Jetz 2012 method).
  dated: true
  dating_method: "Inherited from Pyron 2013 Bayesian dating."
  fossil_calibrations: null
  support_metric: posterior_probability
  divergence_ci_available: true
  divergence_ci_note: "10 posterior tree sets (1000 trees each, ~456 MB total) available on Dryad."
  posterior_n: 10000
  methods_url: https://www.sciencedirect.com/science/article/abs/pii/S0006320716301197

shipped_subclades: []
candidate_subclades:
  - taxon: Iguania (iguanas, anoles, agamids)
    rationale: "~1,800 species; well-studied. Daza et al. 2012 and others have dedicated multi-locus trees."
    candidate_papers: []
    download_target: unknown
    recipe_file: null
    priority: low
    status: proposed
    notes: "Tonini 2016's imputation handles iguanian diversity well; dedicated tree would be supplementary."

resolutions:
  - issue: "Tonini 2016 reports ~9,500 squamate species; data_estimates.csv records 11,744 from Uetz 2025."
    resolved_by: "Both correct for their dates. Tonini 2016 used 2015-era Reptile Database. The 9,755 tip count reflects the tree's vintage; the 11,744 described count reflects Uetz 2025. The 21% gap is real undescribed-since-2016 + cryptic recognition."

website:
  surface_methods: true
  surface_uncertainty: true
  surface_candidates: true
  badge: verified

api_calls:
  - service: crossref
    target: 10.1016/j.biocon.2016.03.039
    ts: 2026-05-18T00:00Z
    outcome: hit
    cache: pending
---

# Squamates — partition audit

**Status: verified.** Tonini 2016 (tree) verified at Crossref. Uetz 2025
(Reptile Database) is a live database — no DOI but a stable URL.

Interval is partly-heuristic: low (11,744) = described count; total (12,000)
and high (13,000) are round-number heuristic placeholders. No paper publishes
a true-diversity range for squamates with the rigour Stork 2018 has for
insects. The tree's 9,755-tip vintage is 2016; described count is 2025.
