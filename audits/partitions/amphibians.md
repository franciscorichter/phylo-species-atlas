---
partition: Amphibians
category: Vertebrates
audit_version: 1
status: verified
last_audited: 2026-05-18
auditor: francisco

estimate_source:
  key: Mora 2011
  doi: 10.1371/journal.pbio.1001127
  doi_verified: true
  crossref_match: true
  year_published: 2011
  cited_value: 15000
  newer_authoritative_known: true
  newer_candidate:
    doi: null
    year: 2025
    study: AmphibiaWeb live database
    cited_value: 8863
    recommendation: watch
    notes: "AmphibiaWeb tracks described species count continuously (~8,863 in 2025); Mora 2011 estimated true diversity at 15,000."

canonical_tree:
  filename: amphibians.nwk
  provenance_group: amphibians
  study: Jetz W, Pyron RA (2018)
  doi: 10.1038/s41559-018-0515-5
  doi_verified: true
  crossref_match: true
  newer_authoritative_known: false
  newer_candidate: null

methods:
  inference: bayesian
  software: ["BEAST (Pyron 2014 base)", "Jetz 2012 imputation"]
  data_type: supermatrix
  n_loci: 15
  dated: true
  dating_method: "15-gene supermatrix Bayesian + taxonomic imputation. Covers frogs + salamanders + caecilians."
  support_metric: posterior_probability
  divergence_ci_available: true
  divergence_ci_note: "10 posterior tree sets available on Dryad."
  posterior_n: null

shipped_subclades:
  - filename: frogs.nwk
    study: Portik DM, Streicher JW, Wiens JJ (2023)
    taxon: Anura (frogs)
    ntips: 5326
    methods_inherits_parent: false
  - filename: salamanders.nwk
    study: Stewart AA, Wiens JJ (2025)
    taxon: Caudata (salamanders)
    ntips: 796
    methods_inherits_parent: false

candidate_subclades:
  - taxon: Gymnophiona (caecilians)
    rationale: "~215 caecilian species. Both Jetz-Pyron 2018 and the shipped sub-clades cover Anura + Caudata. A dedicated caecilian tree (San Mauro 2014) would complete the order-level coverage."
    candidate_papers:
      - doi: null
        year: 2014
        study: San Mauro et al. — caecilian phylogeny
        ntips: null
        url: null
    download_target: unknown
    recipe_file: null
    priority: low
    status: proposed

resolutions:
  - issue: "Mora 2011 estimated 15,000 species; current AmphibiaWeb count is 8,863. Gap reflects undescribed cryptic diversity, not a counting error."
    resolved_by: "Documented. Mora 2011's estimate is bilateral (10,000-20,000); paper-published range. Honest scientific uncertainty."

website:
  surface_methods: true
  surface_uncertainty: true
  surface_candidates: true
  badge: verified
---

# Amphibians — partition audit

**Status: verified.** Jetz-Pyron 2018 tree paper DOI-verified. Mora 2011
publishes a true paper-derived bilateral estimate (10,000-20,000), unlike
checklist-based partitions where the high bound is heuristic. Coverage:
7,239 / 8,863 = 81.7% of described, but only ~48% of Mora 2011's central
estimate — substantial undescribed amphibian diversity.

Sub-clades: frogs (Portik 2023, 5,326 tips) + salamanders (Stewart 2025,
796 tips) shipped. Caecilians (~215 species) under-represented; a dedicated
San Mauro 2014-style tree would complete order-level coverage.
