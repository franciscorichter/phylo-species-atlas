---
partition: Spiders
category: Arthropods
audit_version: 1
status: verified
last_audited: 2026-05-18
auditor: francisco

estimate_source:
  key: Agnarsson 2013
  doi: null
  doi_verified: unchecked
  year_published: 2013
  cited_value: 120000

canonical_tree:
  filename: spiders.nwk
  provenance_group: spiders
  study: Garrison NL et al. (2016)
  doi: 10.7717/peerj.1719
  doi_verified: true
  crossref_match: true

methods:
  inference: ml
  software: [RAxML]
  data_type: UCE
  dating_method: "UCEs + transcriptomes, time-calibrated."
  support_metric: bootstrap
  divergence_ci_available: true
  note: "Alternative calibration and ML backbone also available."

shipped_subclades: []
candidate_subclades: []

resolutions:
  - issue: "1,456-tip tree covers 2.8% of described spider species (~52k). World Spider Catalog 2025 tracks the canonical described count continuously."
    resolved_by: "Documented. Sparse coverage reflects state of spider phylogenetics; no comprehensive tree exists."

website:
  surface_methods: true
  surface_uncertainty: true
  surface_candidates: true
  badge: verified
---

# Spiders — partition audit

**Status: verified.** Garrison 2016 PeerJ (open access) tree (DOI-verified)
covers 1,456 species via UCEs + transcriptomes. Agnarsson 2013 estimate
(80k-170k range) appears to be a true paper-published bilateral interval.
World Spider Catalog tracks live described count.
