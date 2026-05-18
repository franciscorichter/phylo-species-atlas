---
partition: Fish
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
  cited_value: 40000
  newer_authoritative_known: true
  newer_candidate:
    doi: null
    year: 2025
    study: Eschmeyer's Catalog of Fishes (live)
    cited_value: 36000
    recommendation: watch
    notes: "Eschmeyer's Catalog tracks the canonical described count of ray-finned fishes."

canonical_tree:
  filename: fish.nwk
  provenance_group: fish
  study: Rabosky DL, Chang J, Title PO, Cowman PF, Sallan L, Friedman M et al. (2018)
  doi: 10.1038/s41586-018-0273-1
  doi_verified: true
  crossref_match: true
  newer_authoritative_known: false
  newer_candidate: null

methods:
  inference: ml
  software: [RAxML, treePL]
  data_type: supermatrix
  n_loci: 27
  dated: true
  dating_method: "27-gene supermatrix RAxML ML, treePL dating with 130 fossil calibrations."
  fossil_calibrations: 130
  support_metric: bootstrap
  divergence_ci_available: true
  divergence_ci_note: "11,638-tip dated backbone + 100 chronograms with taxonomic imputation (FishTreeOfLife)."
  posterior_n: 100
  methods_url: https://fishtreeoflife.org/

shipped_subclades:
  - filename: neotropical_fish.nwk
    study: Tagliacollo VA et al. (2024)
    taxon: Neotropical freshwater fishes
    ntips: 3167
    methods_inherits_parent: false

candidate_subclades:
  - taxon: Cichlidae (cichlids)
    rationale: "~1,700 species; iconic radiation. Multiple cichlid trees available; would complement the broad Rabosky tree."
    candidate_papers: []
    priority: low
    status: proposed

resolutions:
  - issue: "Atlas ships the 31,516-tip FishTreeOfLife chronogram with taxonomic imputation; only 11,638 tips have direct DNA data."
    resolved_by: "Documented in methods.note. Users interpreting species-level placement uncertainty should consult Rabosky 2018 §Methods."

website:
  surface_methods: true
  surface_uncertainty: true
  surface_candidates: true
  badge: verified
---

# Fish — partition audit

**Status: verified.** Rabosky 2018 (Nature, DOI-verified) is the canonical
ray-finned fish phylogeny. Mora 2011 publishes a bilateral diversity range
(36,000-50,000) — another paper-published-range case (rare among
vertebrates).

Coverage is fish-tree-of-life vintage: 11,638 tips with direct DNA + ~20k
added via taxonomic imputation. The 31,516-tip atlas tree reflects the
fully-imputed chronogram. Sub-clade: neotropical_fish (Tagliacollo 2024,
3,167 tips).
